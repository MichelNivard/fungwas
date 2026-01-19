// stage1_cpp.cpp - pybind11 wrapper for stage1 block scores
// Compile with: g++ -O3 -Wall -shared -std=c++17 -fPIC $(python3 -m pybind11 --includes) stage1_cpp.cpp -o stage1_cpp$(python3-config --extension-suffix) -larmadillo -fopenmp

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <armadillo>
#include <omp.h>
#include <cmath>

namespace py = pybind11;

// Returns a 3D array: (M, n_blocks, 1+T) 
// If SNP is invalid, all values will be NaN
py::array_t<double> stage1_block_scores(
    py::array_t<double> G_arr,
    py::array_t<double> Y_resid_arr,
    py::array_t<double> Q_arr,
    py::array_t<int> block_ids_arr,
    int n_blocks,
    double min_mac = 5.0,
    double min_var = 1e-8,
    int n_threads = 8
) {
    // Get buffer info
    auto G_buf = G_arr.request();
    auto Y_buf = Y_resid_arr.request();
    auto Q_buf = Q_arr.request();
    auto block_buf = block_ids_arr.request();
    
    // Dimensions
    size_t N = G_buf.shape[0];
    size_t M = G_buf.shape[1];
    size_t T = Y_buf.shape[1];
    size_t K = Q_buf.shape[1];
    
    // Create armadillo views (no copy)
    double* G_ptr = (double*)G_buf.ptr;
    double* Y_ptr = (double*)Y_buf.ptr;
    double* Q_ptr = (double*)Q_buf.ptr;
    int* block_ids = (int*)block_buf.ptr;
    
    // Pre-allocate output: M x n_blocks x (1+T)
    size_t out_cols = 1 + T;
    py::array_t<double> result({(int)M, n_blocks, (int)out_cols});
    auto result_buf = result.request();
    double* out_ptr = (double*)result_buf.ptr;
    
    // Initialize with NaN
    size_t total_size = M * n_blocks * out_cols;
    for (size_t i = 0; i < total_size; i++) {
        out_ptr[i] = std::nan("");
    }
    
    #ifdef _OPENMP
    omp_set_num_threads(n_threads);
    #endif
    
    // Process each SNP in parallel
    #pragma omp parallel for schedule(dynamic, 20)
    for (size_t j = 0; j < M; j++) {
        // Work with raw pointers for speed
        double g_sum = 0.0;
        size_t n_valid = 0;
        
        // Thread-local storage
        std::vector<double> g(N);
        
        for (size_t i = 0; i < N; i++) {
            double val = G_ptr[i + j * N];  // Column-major access
            g[i] = val;
            if (std::isfinite(val)) {
                g_sum += val;
                n_valid++;
            }
        }
        
        if (n_valid < 100) continue;  // Leave as NaN
        
        // Mean impute
        double g_mean = g_sum / n_valid;
        for (size_t i = 0; i < N; i++) {
            if (!std::isfinite(g[i])) g[i] = g_mean;
        }
        
        // Residualize: g_resid = g - Q * (Q' * g)
        std::vector<double> Qtg(K, 0.0);
        for (size_t k = 0; k < K; k++) {
            double sum = 0.0;
            for (size_t i = 0; i < N; i++) {
                sum += Q_ptr[i + k * N] * g[i];
            }
            Qtg[k] = sum;
        }
        
        std::vector<double> g_resid(N);
        for (size_t i = 0; i < N; i++) {
            double proj = 0.0;
            for (size_t k = 0; k < K; k++) {
                proj += Q_ptr[i + k * N] * Qtg[k];
            }
            g_resid[i] = g[i] - proj;
        }
        
        // Accumulate block stats
        std::vector<double> stats(n_blocks * out_cols, 0.0);
        
        for (size_t i = 0; i < N; i++) {
            int b = block_ids[i];
            if (b < 0 || b >= n_blocks) continue;
            
            double g_i = g_resid[i];
            stats[b * out_cols + 0] += g_i * g_i;
            for (size_t t = 0; t < T; t++) {
                stats[b * out_cols + 1 + t] += g_i * Y_ptr[i + t * N];
            }
        }
        
        // Copy to output
        size_t out_offset = j * n_blocks * out_cols;
        for (size_t bi = 0; bi < (size_t)n_blocks; bi++) {
            for (size_t c = 0; c < out_cols; c++) {
                out_ptr[out_offset + bi * out_cols + c] = stats[bi * out_cols + c];
            }
        }
    }
    
    return result;
}

PYBIND11_MODULE(stage1_cpp, m) {
    m.doc() = "C++ Stage1 block scores with pybind11";
    m.def("stage1_block_scores", &stage1_block_scores,
          "Compute block-level sufficient statistics for jackknife SE estimation.\n"
          "Returns 3D array (M, n_blocks, 1+T). NaN for invalid SNPs.",
          py::arg("G"),
          py::arg("Y_resid"),
          py::arg("Q"),
          py::arg("block_ids"),
          py::arg("n_blocks"),
          py::arg("min_mac") = 5.0,
          py::arg("min_var") = 1e-8,
          py::arg("n_threads") = 8
    );
}
