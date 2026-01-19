// Stage 1 C++ kernel with pybind11 bindings
// High-performance block score computation for jackknife SE estimation

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <armadillo>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace py = pybind11;

// Compute block scores for jackknife covariance estimation
// Returns: (M, n_blocks, 1+T) array of per-block statistics
py::array_t<double> stage1_block_scores(
    py::array_t<double, py::array::f_style | py::array::forcecast> G_arr,
    py::array_t<double, py::array::f_style | py::array::forcecast> RIF_resid_arr,
    py::array_t<double, py::array::f_style | py::array::forcecast> Q_arr,
    py::array_t<int32_t> block_ids_arr,
    int n_blocks,
    double min_mac,
    double min_var,
    int n_threads
) {
    // Get array info
    auto G_buf = G_arr.request();
    auto RIF_buf = RIF_resid_arr.request();
    auto Q_buf = Q_arr.request();
    auto block_buf = block_ids_arr.request();
    
    size_t N = G_buf.shape[0];
    size_t M = G_buf.shape[1];
    size_t T = RIF_buf.shape[1];
    size_t K = Q_buf.shape[1];
    
    double* G_ptr = static_cast<double*>(G_buf.ptr);
    double* RIF_ptr = static_cast<double*>(RIF_buf.ptr);
    double* Q_ptr = static_cast<double*>(Q_buf.ptr);
    int32_t* block_ptr = static_cast<int32_t*>(block_buf.ptr);
    
    // Wrap in Armadillo (no copy for column-major)
    arma::mat G(G_ptr, N, M, false, true);
    arma::mat RIF_resid(RIF_ptr, N, T, false, true);
    arma::mat Q(Q_ptr, N, K, false, true);
    
    // Output array: (M, n_blocks, 1+T)
    std::vector<size_t> out_shape = {M, (size_t)n_blocks, 1 + T};
    py::array_t<double> result(out_shape);
    auto result_buf = result.request();
    double* result_ptr = static_cast<double*>(result_buf.ptr);
    
    // Initialize to NaN
    size_t total = M * n_blocks * (1 + T);
    for (size_t i = 0; i < total; i++) {
        result_ptr[i] = std::numeric_limits<double>::quiet_NaN();
    }
    
    #ifdef _OPENMP
    omp_set_num_threads(n_threads);
    #endif
    
    #pragma omp parallel for schedule(dynamic, 20)
    for (size_t j = 0; j < M; j++) {
        arma::vec g = G.col(j);
        
        // Count valid and compute mean
        size_t n_valid = 0;
        double g_sum = 0.0;
        for (size_t i = 0; i < N; i++) {
            if (std::isfinite(g(i))) {
                g_sum += g(i);
                n_valid++;
            }
        }
        
        if (n_valid < 100) continue;
        
        double g_mean = g_sum / n_valid;
        double mac = std::min(g_sum, 2.0 * n_valid - g_sum);
        if (mac < min_mac) continue;
        
        // Mean impute missing
        for (size_t i = 0; i < N; i++) {
            if (!std::isfinite(g(i))) g(i) = g_mean;
        }
        
        // Residualize on covariates
        arma::vec Qtg = Q.t() * g;
        arma::vec g_resid = g - Q * Qtg;
        
        double g2 = arma::dot(g_resid, g_resid);
        if (g2 < min_var) continue;
        
        // Accumulate block scores
        arma::mat stats(n_blocks, 1 + T, arma::fill::zeros);
        
        for (size_t i = 0; i < N; i++) {
            int b = block_ptr[i];
            if (b < 0 || b >= n_blocks) continue;
            
            double g_i = g_resid(i);
            stats(b, 0) += g_i * g_i;
            for (size_t t = 0; t < T; t++) {
                stats(b, 1 + t) += g_i * RIF_resid(i, t);
            }
        }
        
        // Copy to output
        for (int b = 0; b < n_blocks; b++) {
            for (size_t t = 0; t < 1 + T; t++) {
                // result[j, b, t] in row-major-ish indexing
                size_t idx = j * n_blocks * (1 + T) + b * (1 + T) + t;
                result_ptr[idx] = stats(b, t);
            }
        }
    }
    
    return result;
}

// Simple info function
py::dict openmp_info() {
    py::dict info;
    #ifdef _OPENMP
    info["openmp_enabled"] = true;
    info["max_threads"] = omp_get_max_threads();
    #else
    info["openmp_enabled"] = false;
    info["max_threads"] = 1;
    #endif
    return info;
}

PYBIND11_MODULE(_stage1_cpp, m) {
    m.doc() = "FungWas Stage 1 C++ kernel";
    
    m.def("stage1_block_scores", &stage1_block_scores,
          "Compute block scores for jackknife covariance",
          py::arg("G"), py::arg("RIF_resid"), py::arg("Q"),
          py::arg("block_ids"), py::arg("n_blocks"),
          py::arg("min_mac") = 5.0, py::arg("min_var") = 1e-8,
          py::arg("n_threads") = 1);
    
    m.def("openmp_info", &openmp_info, "Get OpenMP configuration");
}
