#ifdef TIME

#include <sstream>
#include <vector>

#include <omp.h>

#include "common.h"
#include "io_enhanced.h"
#include "timer.h"

int timer_monitoring_state;

std::vector<std::vector<double> > compute_step_kernel_times;
std::vector<std::vector<double> > compute_flux_edge_kernel_times;
std::vector<std::vector<double> > update_kernel_times;
std::vector<std::vector<double> > time_step_kernel_times;
std::vector<std::vector<double> > restrict_kernel_times;
std::vector<std::vector<double> > prolong_kernel_times;
std::vector<std::vector<double> > unstructured_stream_kernel_times;
std::vector<std::vector<double> > unstructured_compute_kernel_times;

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

    compute_step_kernel_times.resize(num_threads);
    compute_flux_edge_kernel_times.resize(num_threads);
    update_kernel_times.resize(num_threads);
    time_step_kernel_times.resize(num_threads);
    restrict_kernel_times.resize(num_threads);
    prolong_kernel_times.resize(num_threads);
    unstructured_stream_kernel_times.resize(num_threads);
    unstructured_compute_kernel_times.resize(num_threads);

    start_times.resize(num_threads);

    for (int t=0; t<num_threads; t++) {
        compute_step_kernel_times.at(t).resize(levels, 0.0f);
        compute_flux_edge_kernel_times.at(t).resize(levels, 0.0f);
        update_kernel_times.at(t).resize(levels, 0.0f);
        time_step_kernel_times.at(t).resize(levels, 0.0f);
        restrict_kernel_times.at(t).resize(levels, 0.0f);
        prolong_kernel_times.at(t).resize(levels, 0.0f);
        unstructured_stream_kernel_times.at(t).resize(levels, 0.0f);
        unstructured_compute_kernel_times.at(t).resize(levels, 0.0f);

        start_times.at(t) = 0.0f;
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

    if (current_kernel == COMPUTE_STEP) {
        compute_step_kernel_times[tid][level] += duration;
    }
    else if (current_kernel == COMPUTE_FLUX_EDGE) {
        compute_flux_edge_kernel_times[tid][level] += duration;
    }
    else if (current_kernel == UPDATE) {
        update_kernel_times[tid][level] += duration;
    }
    else if (current_kernel == TIME_STEP) {
        time_step_kernel_times[tid][level] += duration;
    }
    else if (current_kernel == RESTRICT) {
        restrict_kernel_times[tid][level] += duration;
    }
    else if (current_kernel == PROLONG) {
        prolong_kernel_times[tid][level] += duration;
    }
    else if (current_kernel == UNSTRUCTURED_STREAM) {
        unstructured_stream_kernel_times[tid][level] += duration;
    }
    else if (current_kernel == UNSTRUCTURED_COMPUTE) {
        unstructured_compute_kernel_times[tid][level] += duration;
    }
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
        data_line << compute_flux_edge_kernel_times[tid][l] << "," ;
        data_line << update_kernel_times[tid][l] << "," ;
        data_line << compute_step_kernel_times[tid][l] << "," ;
        data_line << time_step_kernel_times[tid][l] << "," ;
        data_line << restrict_kernel_times[tid][l] << "," ;
        data_line << prolong_kernel_times[tid][l] << "," ;
        data_line << unstructured_stream_kernel_times[tid][l] << "," ;
        data_line << unstructured_compute_kernel_times[tid][l] << "," ;
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
