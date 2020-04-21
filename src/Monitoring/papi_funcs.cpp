#ifdef PAPI

#include <sstream>
#include <string.h>
#include <vector>

#include "common.h"
#include "io_enhanced.h"
#include "papi_funcs.h"

int papi_monitoring_state;

std::vector<int> event_sets;
std::vector<int> num_thread_events;
std::vector<std::vector<int> > thread_events;

std::vector<std::vector<long_long> > compute_step_kernel_event_counts;
std::vector<std::vector<long_long> > compute_flux_edge_kernel_event_counts;
std::vector<std::vector<long_long> > update_kernel_event_counts;
std::vector<std::vector<long_long> > indirect_rw_kernel_event_counts;
std::vector<std::vector<long_long> > time_step_kernel_event_counts;
std::vector<std::vector<long_long> > restrict_kernel_event_counts;
std::vector<std::vector<long_long> > prolong_kernel_event_counts;

std::vector<long_long*> temp_count_stores;

void init_papi()
{
    int ret;

    log("Initialising PAPI");

    papi_monitoring_state = 1;

    #ifdef OMP
    int num_threads = omp_get_max_threads();
    #else
    int num_threads = 1;
    #endif

    num_thread_events.resize(num_threads, 0);
    event_sets.resize(num_threads, PAPI_NULL);
    thread_events.resize(num_threads);

    compute_step_kernel_event_counts.resize(num_threads);
    compute_flux_edge_kernel_event_counts.resize(num_threads);
    update_kernel_event_counts.resize(num_threads);
    indirect_rw_kernel_event_counts.resize(num_threads);
    time_step_kernel_event_counts.resize(num_threads);
    restrict_kernel_event_counts.resize(num_threads);
    prolong_kernel_event_counts.resize(num_threads);

    temp_count_stores.resize(num_threads);

    if (PAPI_num_counters() < 2) {
       fprintf(stderr, "No hardware counters here, or PAPI not supported.\n");
       exit(-1);
    }
    log("PAPI_num_counters() complete");

    if ((ret=PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) {
        fprintf(stderr, "PAPI_library_init() failed: '%s'.\n", PAPI_strerror(ret));
        exit(EXIT_FAILURE);
    }
    log("PAPI_library_init() complete");

    #ifdef OMP
        #pragma omp parallel
        {
            #pragma omp critical
            {
                // if ((ret=PAPI_thread_init(pthread_self)) != PAPI_OK) {
                if ((ret=PAPI_thread_init(omp_get_thread_num_ul)) != PAPI_OK) {
                    fprintf(stderr, "Thread %d: PAPI_thread_init() failed: '%s'.\n", omp_get_thread_num(), PAPI_strerror(ret));
                    exit(EXIT_FAILURE);
                }
            }
        }
    #else
        // if ((ret=PAPI_thread_init(pthread_self)) != PAPI_OK) {
        if ((ret=PAPI_thread_init(omp_get_thread_num_ul)) != PAPI_OK) {
            fprintf(stderr, "PAPI_thread_init() failed: '%s'.\n", PAPI_strerror(ret));
            exit(EXIT_FAILURE);
        }
    #endif
    log("PAPI_thread_init() complete");
}

void load_papi_events()
{
    int ret;

    bool first_thread_is_initialising = true;
    #if defined OMP
        #pragma omp parallel private(ret)
        {
            const int tid = omp_get_thread_num();
            #pragma omp critical
            {
    #else
        const int tid = 0;
    #endif
        event_sets[tid] = PAPI_NULL;
        ret = PAPI_create_eventset(&event_sets[tid]);
        if (ret != PAPI_OK || event_sets[tid]==PAPI_NULL) {
            fprintf(stderr, "Thread %d on CPU %d: PAPI_create_eventset() failed: '%s'.\n", tid, sched_getcpu(), PAPI_strerror(ret));
            DEBUGGABLE_ABORT
        }

        char event_name[512];
        std::string line;
        std::ifstream file_if(conf.papi_config_file);
        if(!file_if.is_open()) {
            printf("ERROR: Failed to open PAPI config file: '%s'\n", conf.papi_config_file);
            exit(EXIT_FAILURE);
        }
        while(std::getline(file_if, line))
        {
            if (line.c_str()[0] == '#' || strcmp(line.c_str(), "")==0) {
                continue;
            }

            if (line.find(std::string("unc_imc")) != std::string::npos) {
                // This is a uncore imc event. If many threads monitor it 
                // then PAPI throws up errors, so allow only one thread to 
                // monitor it.
                if (!first_thread_is_initialising) {
                    continue;
                }
            }
            if (line.find(std::string("unc_edc")) != std::string::npos) {
                // This is a uncore mcdram event. If many threads monitor it 
                // then PAPI throws up errors, so allow only one thread to 
                // monitor it.
                if (!first_thread_is_initialising) {
                  continue;
                }
            }

            strcpy(event_name, line.c_str());

            ret = PAPI_add_named_event(event_sets[tid], event_name);
            // ret = PAPI_OK - 1;
            if (ret != PAPI_OK) {
                // printf("PAPI_add_named_event() failed, attempting add by code conversion\n");

                // fprintf(stderr, "Thread %d on CPU %d: failed to add event %s to event set: '%s'.\n", omp_get_thread_num(), sched_getcpu(), event_name, PAPI_strerror(ret));
                // if (event_sets[tid]==PAPI_NULL) {
                //     fprintf(stderr, "... and event_set=PAPI_NULL\n");
                // }

                // It appears that PAPI_add_named_event() only works with native events, not PAPI presets.
                int code = -1;
                ret = PAPI_event_name_to_code(event_name, &code);
                if (ret != PAPI_OK) {
                    printf("Could not convert string '%s' to PAPI event, error = %s\n", event_name, PAPI_strerror(ret));
                } else {
                    if (PAPI_query_event(code) != PAPI_OK) {
                        printf("PAPI event %s not present\n", event_name);
                    } else {
                        ret = PAPI_add_event(event_sets[tid], code);
                        if (ret != PAPI_OK) {
                            #if defined OMP
                                fprintf(stderr, "Thread %d on CPU %d: failed to add event %d to event set: '%s'.\n", omp_get_thread_num(), sched_getcpu(), code, PAPI_strerror(ret));
                            #else
                                fprintf(stderr, "Failed to add event %d to event set: '%s'.\n", code, PAPI_strerror(ret));
                            #endif
                            if (event_sets[tid]==PAPI_NULL) {
                                fprintf(stderr, "... and event_set=PAPI_NULL\n");
                            }
                            DEBUGGABLE_ABORT
                        }
                        else {
                            log("Monitoring PAPI event '%s'\n", event_name);
                        }
                    }
                }
            }
        }
        if (file_if.bad()) {
            printf("ERROR: Failed to read PAPI config file: %s\n", conf.papi_config_file);
            exit(EXIT_FAILURE);
        }
        log("Finished parsing PAPI file");

        int n = PAPI_num_events(event_sets[tid]);
        num_thread_events[tid] = n;
        if (n == 0) {
            event_sets[tid] = PAPI_NULL;
        }
        
        first_thread_is_initialising = false;

    #if defined OMP
            } // End of critical region
        } // End of parallel loop
    #endif

    #ifdef OMP
    int num_threads = omp_get_max_threads();
    #else
    int num_threads = 1;
    #endif
    for (int tid=0; tid<num_threads; tid++) {
        int n = num_thread_events[tid];

        compute_step_kernel_event_counts.at(tid).resize(n*levels);
        compute_flux_edge_kernel_event_counts.at(tid).resize(n*levels);
        update_kernel_event_counts.at(tid).resize(n*levels);
        indirect_rw_kernel_event_counts.at(tid).resize(n*levels);
        time_step_kernel_event_counts.at(tid).resize(n*levels);
        restrict_kernel_event_counts.at(tid).resize(n*levels);
        prolong_kernel_event_counts.at(tid).resize(n*levels);

        temp_count_stores.at(tid) = alloc<long_long>(n*levels);

        thread_events.at(tid).resize(n);
        int* temp_thread_events = alloc<int>(n);
        if (PAPI_list_events(event_sets.at(tid), temp_thread_events, &n) != PAPI_OK) {
            fprintf(stderr, "ERROR: PAPI_list_events() failed\n");
            DEBUGGABLE_ABORT
        }
        for (int i=0; i<n; i++) {
            thread_events.at(tid).at(i) = temp_thread_events[i];
        }
        dealloc(temp_thread_events);
    }

    #if defined OMP
        #pragma omp parallel
        {
        const int tid = omp_get_thread_num();
    #endif
        if (event_sets[tid] != PAPI_NULL) {
            if (PAPI_start(event_sets[tid]) != PAPI_OK) {
                fprintf(stderr, "ERROR: Thread %d failed to start PAPI\n", tid);
                DEBUGGABLE_ABORT
            }
            if (PAPI_stop(event_sets[tid], temp_count_stores[tid]) != PAPI_OK) {
                fprintf(stderr, "ERROR: Thread %d failed to stop PAPI\n", tid);
                DEBUGGABLE_ABORT
            }
        }
    #if defined OMP
        } // End of parallel loop
    #endif
}


void start_papi()
{
    if (!papi_monitoring_state) return;

    #ifdef OMP
        const int tid = omp_get_thread_num();
    #else
        const int tid = 0;
    #endif

    if (PAPI_start(event_sets[tid]) != PAPI_OK) {
        fprintf(stderr, "ERROR: Thread %d failed to start PAPI\n", tid);
        exit(EXIT_FAILURE);
    }
}

void stop_papi()
{
    if (!papi_monitoring_state) return;

    #ifdef OMP
        const int tid = omp_get_thread_num();
    #else
        const int tid = 0;
    #endif

    if (PAPI_stop(event_sets[tid], temp_count_stores[tid]) != PAPI_OK) {
        fprintf(stderr, "ERROR: Thread %d failed to stop PAPI\n", tid);
        exit(EXIT_FAILURE);
    }

    int num_events = num_thread_events[tid];
    if (current_kernel == COMPUTE_STEP) {
        for (int e=0; e<num_events; e++) {
            compute_step_kernel_event_counts[tid][level*num_events + e] += temp_count_stores[tid][e];
        }
    }
    else if (current_kernel == COMPUTE_FLUX_EDGE) {
        for (int e=0; e<num_events; e++) {
            compute_flux_edge_kernel_event_counts[tid][level*num_events + e] += temp_count_stores[tid][e];
        }
    }
    else if (current_kernel == UPDATE) {
        for (int e=0; e<num_events; e++) {
            update_kernel_event_counts[tid][level*num_events + e] += temp_count_stores[tid][e];
        }
    }
    else if (current_kernel == INDIRECT_RW) {
        for (int e=0; e<num_events; e++) {
            indirect_rw_kernel_event_counts[tid][level*num_events + e] += temp_count_stores[tid][e];
        }
    }
    else if (current_kernel == TIME_STEP) {
        for (int e=0; e<num_events; e++) {
            time_step_kernel_event_counts[tid][level*num_events + e] += temp_count_stores[tid][e];
        }
    }
    else if (current_kernel == RESTRICT) {
        for (int e=0; e<num_events; e++) {
            restrict_kernel_event_counts[tid][level*num_events + e] += temp_count_stores[tid][e];
        }
    }
    else if (current_kernel == PROLONG) {
        for (int e=0; e<num_events; e++) {
            prolong_kernel_event_counts[tid][level*num_events + e] += temp_count_stores[tid][e];
        }
    }
}


void dump_papi_counters_to_file(int size)
{
    log("Called dump_papi_counters_to_file()");

    // File admin:
    std::string filepath = std::string(conf.output_file_prefix);
    if (filepath.size() > 0 && filepath.at(filepath.size()-1) != '/') {
        filepath += ".";
    }
    filepath += "PAPI.csv";

    if (file_exists(filepath.c_str())) {
        std::remove(filepath.c_str());
    }

    bool write_header = !file_exists(filepath.c_str());

    std::ofstream outfile;
    if (write_header) {
        outfile.open(filepath.c_str(), std::ios_base::out);
    } else {
        outfile.open(filepath.c_str(), std::ios_base::app);
    }

    std::string header_s;
    std::string data_line_s;
    prepare_csv_identification(write_header, &header_s, &data_line_s, size);

    if (write_header) {
        std::ostringstream header;
        header << header_s;
        header << "ThreadNum,";
        header << "CpuId,";
        header << "PAPI counter,";
        for (int l=0; l<levels; l++) {
            header << "flux" << l << "," ;
            header << "update" << l << "," ;
            header << "compute_step" << l << "," ;
            header << "time_step" << l << "," ;
            header << "restrict" << l << "," ;
            header << "prolong" << l << "," ;
            header << "indirect_rw" << l << "," ;
        }
        outfile << header.str() << std::endl;
    }

    #ifdef OMP
        int cpu_ids[conf.omp_num_threads];
        #pragma omp parallel
        {
            cpu_ids[omp_get_thread_num()] = sched_getcpu();
        }
        for (int tid=0; tid<conf.omp_num_threads; tid++) {
            const int cpu_id = cpu_ids[tid];
    #else
        const int cpu_id = sched_getcpu();
        int tid = 0;
    #endif

    int num_events = num_thread_events[tid];

    for (int eid=0; eid<num_thread_events[tid]; eid++)
    {
        std::ostringstream event_data_line;
        event_data_line << data_line_s;

        char eventName[PAPI_MAX_STR_LEN] = "";
        if (PAPI_event_code_to_name(thread_events[tid][eid], eventName) != PAPI_OK) {
            fprintf(stderr, "ERROR: Failed to convert code %d to name (tid=%d)\n", thread_events[tid][eid], tid);
            DEBUGGABLE_ABORT
        }

        event_data_line << tid << ",";
        event_data_line << cpu_id << ",";

        event_data_line << eventName << ",";

        for (int l=0; l<levels; l++) {
            const int idx = l*num_events + eid;
            event_data_line << compute_flux_edge_kernel_event_counts[tid][idx] << ',';
            event_data_line << update_kernel_event_counts[tid][idx] << ',';
            event_data_line << compute_step_kernel_event_counts[tid][idx] << ',';
            event_data_line << time_step_kernel_event_counts[tid][idx] << ',';
            event_data_line << restrict_kernel_event_counts[tid][idx] << ',';
            event_data_line << prolong_kernel_event_counts[tid][idx] << ',';
            event_data_line << indirect_rw_kernel_event_counts[tid][idx] << ',';
        }
        
        outfile << event_data_line.str() << std::endl;
    }

    #if defined OMP
        } // End loop over thread data.
    #endif

    outfile.close();

    if (!file_exists(filepath.c_str())) {
        printf("Failed to write PAPI counts to: %s\n", filepath.c_str());
    } else {
        printf("PAPI counts written to: %s\n", filepath.c_str());
    }
}

#endif
