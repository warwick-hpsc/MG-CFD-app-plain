#include <sstream>
#include <vector>

#include "io_enhanced.h"
#include "loop_stats.h"

int iters_monitoring_state;

std::vector<std::vector<long> > compute_step_kernel_niters;
std::vector<std::vector<long> > compute_flux_edge_kernel_niters;
std::vector<std::vector<long> > update_kernel_niters;
std::vector<std::vector<long> > time_step_kernel_niters;
std::vector<std::vector<long> > restrict_kernel_niters;
std::vector<std::vector<long> > prolong_kernel_niters;
std::vector<std::vector<long> > unstructured_stream_kernel_niters;

void init_iters()
{
    log("Initialising iters");

    #ifdef OMP
    int num_threads = omp_get_max_threads();
    #else
    int num_threads = 1;
    #endif

    iters_monitoring_state = 1;

    compute_step_kernel_niters.resize(num_threads);
    compute_flux_edge_kernel_niters.resize(num_threads);
    update_kernel_niters.resize(num_threads);
    time_step_kernel_niters.resize(num_threads);
    restrict_kernel_niters.resize(num_threads);
    prolong_kernel_niters.resize(num_threads);
    unstructured_stream_kernel_niters.resize(num_threads);

    for (int t=0; t<num_threads; t++) {
        compute_step_kernel_niters.at(t).resize(levels, 0);
        compute_flux_edge_kernel_niters.at(t).resize(levels, 0);
        update_kernel_niters.at(t).resize(levels, 0);
        time_step_kernel_niters.at(t).resize(levels, 0);
        restrict_kernel_niters.at(t).resize(levels, 0);
        prolong_kernel_niters.at(t).resize(levels, 0);
        unstructured_stream_kernel_niters.at(t).resize(levels, 0);
    }
}

void record_iters(long loop_start, long loop_end)
{
    if (!iters_monitoring_state) return;

    #ifdef OMP
        const int tid = omp_get_thread_num();
    #else
        const int tid = 0;
    #endif

    long niters = loop_end - loop_start;

    if (current_kernel == COMPUTE_STEP) {
        compute_step_kernel_niters[tid][level] += niters;
    }
    else if (current_kernel == COMPUTE_FLUX_EDGE) {
        compute_flux_edge_kernel_niters[tid][level] += niters;
    }
    else if (current_kernel == UPDATE) {
        update_kernel_niters[tid][level] += niters;
    }
    else if (current_kernel == TIME_STEP) {
        time_step_kernel_niters[tid][level] += niters;
    }
    else if (current_kernel == RESTRICT) {
        restrict_kernel_niters[tid][level] += niters;
    }
    else if (current_kernel == PROLONG) {
        prolong_kernel_niters[tid][level] += niters;
    }
    else if (current_kernel == UNSTRUCTURED_STREAM) {
        unstructured_stream_kernel_niters[tid][level] += niters;
    }
}

void dump_loop_stats_to_file(int size)
{
    log("dump_loop_stats_to_file() called");

    // File admin:
    std::string filepath = std::string(conf.output_file_prefix);
    if (filepath.size() > 0 && filepath.at(filepath.size()-1) != '/') {
        filepath += ".";
    }
    filepath += "LoopNumIters.csv";

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

    std::ostringstream data_line;
    data_line << data_line_s;
    data_line << tid << ",";
    data_line << cpu_id << ",";

    for (int l=0; l<levels; l++) {
        data_line << compute_flux_edge_kernel_niters[tid][l] << "," ;
        data_line << update_kernel_niters[tid][l] << "," ;
        data_line << compute_step_kernel_niters[tid][l] << "," ;
        data_line << time_step_kernel_niters[tid][l] << "," ;
        data_line << restrict_kernel_niters[tid][l] << "," ;
        data_line << prolong_kernel_niters[tid][l] << "," ;
        data_line << unstructured_stream_kernel_niters[tid][l] << "," ;
    }

    outfile << data_line.str() << std::endl;

    #if defined OMP
        } // End loop over thread data.
    #endif

    outfile.close();

    if (!file_exists(filepath.c_str())) {
        printf("Failed to write loop stats to: %s\n", filepath.c_str());
    } else {
        printf("Loop stats written to: %s\n", filepath.c_str());
    }

    log("dump_loop_stats_to_file() complete");
}
