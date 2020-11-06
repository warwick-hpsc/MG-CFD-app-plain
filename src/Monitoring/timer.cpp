#ifdef TIME

#include <sstream>
#include <vector>

#include <omp.h>

#include "common.h"
#include "io_enhanced.h"
#include "timer.h"

int timer_monitoring_state;

std::vector<std::vector<double> > kernel_times[NUM_KERNELS];

std::vector<double> start_times;

void init_timers()
{
    log("Initialising timers");

    timer_monitoring_state = 1;

    #ifdef OMP
    int num_threads = omp_get_max_threads();
    #else
    int num_threads = 1;
    #endif

    start_times.resize(num_threads);
    for (int t=0; t<num_threads; t++) {
        start_times.at(t) = 0.0f;
    }

    for (int nk=0; nk<NUM_KERNELS; nk++) {
        kernel_times[nk].resize(num_threads);
        for (int t=0; t<num_threads; t++) {
            kernel_times[nk].at(t).resize(levels, 0.0f);
        }
    }
}

void start_timer()
{
    if (!timer_monitoring_state) return;

    #ifdef OMP
        const int tid = omp_get_thread_num();
    #else
        const int tid = 0;
    #endif

    start_times[tid] = omp_get_wtime();
}

void stop_timer()
{
    if (!timer_monitoring_state) return;

    #ifdef OMP
        const int tid = omp_get_thread_num();
    #else
        const int tid = 0;
    #endif
    
    double duration = omp_get_wtime() - start_times[tid];

    kernel_times[current_kernel][tid][level] += duration;
}

void dump_timers_to_file(int size, double total_time)
{
    log("dump_timers_to_file() called");

    // File admin:
    std::string filepath = std::string(conf.output_file_prefix);
    if (filepath.size() > 0 && filepath.at(filepath.size()-1) != '/') {
        filepath += ".";
    }
    filepath += "Times.csv";

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
        for (int l=0; l<levels; l++) {
            header << "flux" << l << "," ;
            header << "update" << l << "," ;
            header << "compute_step" << l << "," ;
            header << "time_step" << l << "," ;
            header << "restrict" << l << "," ;
            header << "prolong" << l << "," ;
            header << "unstructured_stream" << l << "," ;
            header << "unstructured_compute" << l << "," ;
            header << "compute_stream" << l << "," ;
        }
        header << "Total," ;
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

    std::ostringstream data_line;
    data_line << data_line_s;
    data_line << tid << ",";
    data_line << cpu_id << ",";

    for (int l=0; l<levels; l++) {
        data_line << kernel_times[COMPUTE_FLUX_EDGE][tid][l] << "," ;
        data_line << kernel_times[UPDATE][tid][l] << "," ;
        data_line << kernel_times[COMPUTE_STEP][tid][l] << "," ;
        data_line << kernel_times[TIME_STEP][tid][l] << "," ;
        data_line << kernel_times[RESTRICT][tid][l] << "," ;
        data_line << kernel_times[PROLONG][tid][l] << "," ;
        data_line << kernel_times[UNSTRUCTURED_STREAM][tid][l] << "," ;
        data_line << kernel_times[UNSTRUCTURED_COMPUTE][tid][l] << "," ;
        data_line << kernel_times[COMPUTE_STREAM][tid][l] << "," ;
    }

    data_line << total_time << "," ;

    outfile << data_line.str() << std::endl;

    #if defined OMP
        } // End loop over thread data.
    #endif

    outfile.close();

    if (!file_exists(filepath.c_str())) {
        printf("Failed to write loop runtimes to: %s\n", filepath.c_str());
    } else {
        printf("Loop runtimes written to: %s\n", filepath.c_str());
    }
}

#endif
