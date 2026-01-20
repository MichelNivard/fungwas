// Stage 1 C++ kernel with pybind11 bindings
// High-performance block score computation for jackknife SE estimation

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <armadillo>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

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

py::list stage1_block_scores_multi(
    py::array_t<double, py::array::f_style | py::array::forcecast> G_arr,
    py::list rif_matrices,
    py::array_t<double, py::array::f_style | py::array::forcecast> Q_arr,
    py::array_t<int32_t> block_ids_arr,
    int n_blocks,
    double min_mac,
    double min_var,
    int n_threads
) {
    auto G_buf = G_arr.request();
    auto Q_buf = Q_arr.request();
    auto block_buf = block_ids_arr.request();

    size_t N = G_buf.shape[0];
    size_t M = G_buf.shape[1];
    size_t K = Q_buf.shape[1];

    if (Q_buf.shape[0] != N) {
        throw std::runtime_error("Q must have same number of rows as G.");
    }

    double* G_ptr = static_cast<double*>(G_buf.ptr);
    double* Q_ptr = static_cast<double*>(Q_buf.ptr);
    int32_t* block_ptr = static_cast<int32_t*>(block_buf.ptr);

    arma::mat G(G_ptr, N, M, false, true);
    arma::mat Q(Q_ptr, N, K, false, true);

    size_t n_pheno = rif_matrices.size();
    py::list results;
    if (n_pheno == 0) {
        return results;
    }

    std::vector<py::array_t<double>> rif_arrays;
    std::vector<arma::mat> rif_mats;
    std::vector<py::array_t<double>> output_arrays;
    std::vector<double*> output_ptrs;
    std::vector<size_t> t_sizes;

    rif_arrays.reserve(n_pheno);
    rif_mats.reserve(n_pheno);
    output_arrays.reserve(n_pheno);
    output_ptrs.reserve(n_pheno);
    t_sizes.reserve(n_pheno);

    double nan_val = std::numeric_limits<double>::quiet_NaN();

    for (size_t p = 0; p < n_pheno; p++) {
        auto rif_arr = rif_matrices[p].cast<py::array_t<double, py::array::f_style | py::array::forcecast>>();
        auto rif_buf = rif_arr.request();
        if (rif_buf.ndim != 2) {
            throw std::runtime_error("Each RIF matrix must be 2D.");
        }
        if (rif_buf.shape[0] != N) {
            throw std::runtime_error("RIF matrices must match G row count.");
        }

        size_t T = rif_buf.shape[1];
        double* rif_ptr = static_cast<double*>(rif_buf.ptr);
        rif_arrays.push_back(rif_arr);
        rif_mats.emplace_back(rif_ptr, N, T, false, true);

        std::vector<size_t> out_shape = {M, (size_t)n_blocks, 1 + T};
        py::array_t<double> out_arr(out_shape);
        auto out_buf = out_arr.request();
        double* out_ptr = static_cast<double*>(out_buf.ptr);
        std::fill(out_ptr, out_ptr + M * n_blocks * (1 + T), nan_val);

        output_arrays.push_back(out_arr);
        output_ptrs.push_back(out_ptr);
        t_sizes.push_back(T);
    }

    #ifdef _OPENMP
    omp_set_num_threads(n_threads);
    #endif

    #pragma omp parallel for schedule(dynamic, 20)
    for (size_t j = 0; j < M; j++) {
        arma::vec g = G.col(j);

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

        for (size_t i = 0; i < N; i++) {
            if (!std::isfinite(g(i))) g(i) = g_mean;
        }

        arma::vec Qtg = Q.t() * g;
        arma::vec g_resid = g - Q * Qtg;

        double g2 = arma::dot(g_resid, g_resid);
        if (g2 < min_var) continue;

        for (size_t p = 0; p < n_pheno; p++) {
            const arma::mat& rif = rif_mats[p];
            size_t T = t_sizes[p];
            arma::mat stats(n_blocks, 1 + T, arma::fill::zeros);

            for (size_t i = 0; i < N; i++) {
                int b = block_ptr[i];
                if (b < 0 || b >= n_blocks) continue;

                double g_i = g_resid(i);
                stats(b, 0) += g_i * g_i;
                for (size_t t = 0; t < T; t++) {
                    stats(b, 1 + t) += g_i * rif(i, t);
                }
            }

            double* out_ptr = output_ptrs[p];
            size_t stride = n_blocks * (1 + T);
            size_t out_offset = j * stride;
            for (int b = 0; b < n_blocks; b++) {
                for (size_t t = 0; t < 1 + T; t++) {
                    size_t idx = out_offset + b * (1 + T) + t;
                    out_ptr[idx] = stats(b, t);
                }
            }
        }
    }

    for (size_t p = 0; p < n_pheno; p++) {
        results.append(output_arrays[p]);
    }

    return results;
}

py::list process_block_dense_impl(
    py::array_t<double, py::array::f_style | py::array::forcecast> G_arr,
    py::list rif_matrices,
    py::array_t<double, py::array::f_style | py::array::forcecast> Q_arr,
    double min_mac,
    double min_var,
    int n_threads
) {
    auto G_buf = G_arr.request();
    auto Q_buf = Q_arr.request();

    size_t N = G_buf.shape[0];
    size_t M = G_buf.shape[1];
    size_t K = Q_buf.shape[1];
    if (Q_buf.shape[0] != N) {
        throw std::runtime_error("Q must have same number of rows as G.");
    }

    double* G_ptr = static_cast<double*>(G_buf.ptr);
    double* Q_ptr = static_cast<double*>(Q_buf.ptr);

    arma::mat G(G_ptr, N, M, false, true);
    arma::mat Q(Q_ptr, N, K, false, true);

    size_t n_pheno = rif_matrices.size();
    py::list results;
    if (n_pheno == 0) {
        return results;
    }

    std::vector<py::array_t<double>> rif_arrays;
    std::vector<arma::mat> rif_mats;
    std::vector<py::array_t<double>> betas_arrays;
    std::vector<py::array_t<double>> se_arrays;
    std::vector<double*> betas_ptrs;
    std::vector<double*> se_ptrs;
    std::vector<size_t> t_sizes;

    rif_arrays.reserve(n_pheno);
    rif_mats.reserve(n_pheno);
    betas_arrays.reserve(n_pheno);
    se_arrays.reserve(n_pheno);
    betas_ptrs.reserve(n_pheno);
    se_ptrs.reserve(n_pheno);
    t_sizes.reserve(n_pheno);

    double nan_val = std::numeric_limits<double>::quiet_NaN();

    for (size_t p = 0; p < n_pheno; p++) {
        auto rif_arr = rif_matrices[p].cast<py::array_t<double, py::array::f_style | py::array::forcecast>>();
        auto rif_buf = rif_arr.request();
        if (rif_buf.ndim != 2) {
            throw std::runtime_error("Each RIF matrix must be 2D.");
        }
        if (rif_buf.shape[0] != N) {
            throw std::runtime_error("RIF matrices must match G row count.");
        }

        size_t T = rif_buf.shape[1];
        double* rif_ptr = static_cast<double*>(rif_buf.ptr);
        rif_arrays.push_back(rif_arr);
        rif_mats.emplace_back(rif_ptr, N, T, false, true);

        py::array_t<double, py::array::f_style> betas_arr({M, T});
        py::array_t<double, py::array::f_style> se_arr({M, T});
        auto betas_buf = betas_arr.request();
        auto se_buf = se_arr.request();
        double* betas_ptr = static_cast<double*>(betas_buf.ptr);
        double* se_ptr = static_cast<double*>(se_buf.ptr);
        std::fill(betas_ptr, betas_ptr + M * T, nan_val);
        std::fill(se_ptr, se_ptr + M * T, nan_val);

        betas_arrays.push_back(betas_arr);
        se_arrays.push_back(se_arr);
        betas_ptrs.push_back(betas_ptr);
        se_ptrs.push_back(se_ptr);
        t_sizes.push_back(T);
    }

    #ifdef _OPENMP
    omp_set_num_threads(n_threads);
    #endif

    #pragma omp parallel for schedule(dynamic, 20)
    for (size_t j = 0; j < M; j++) {
        arma::vec g = G.col(j);

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

        for (size_t i = 0; i < N; i++) {
            if (!std::isfinite(g(i))) g(i) = g_mean;
        }

        arma::vec Qtg = Q.t() * g;
        arma::vec g_resid = g - Q * Qtg;

        double g2 = arma::dot(g_resid, g_resid);
        if (g2 < min_var) continue;

        arma::vec g_sq = g_resid % g_resid;

        for (size_t p = 0; p < n_pheno; p++) {
            const arma::mat& rif = rif_mats[p];
            size_t T = t_sizes[p];
            double* betas_ptr = betas_ptrs[p];
            double* se_ptr = se_ptrs[p];

            for (size_t t = 0; t < T; t++) {
                double g_rif = arma::dot(g_resid, rif.col(t));
                double beta = g_rif / g2;
                betas_ptr[j + M * t] = beta;

                double hc0_sum = 0.0;
                const double* rif_ptr = rif.colptr(t);
                for (size_t i = 0; i < N; i++) {
                    double e_i = rif_ptr[i] - g_resid(i) * beta;
                    hc0_sum += g_sq(i) * e_i * e_i;
                }
                se_ptr[j + M * t] = std::sqrt(hc0_sum) / g2;
            }
        }
    }

    for (size_t p = 0; p < n_pheno; p++) {
        results.append(py::make_tuple(betas_arrays[p], se_arrays[p]));
    }

    return results;
}

py::list process_block_sparse_impl(
    py::array_t<double, py::array::f_style | py::array::forcecast> G_arr,
    py::list rif_matrices,
    py::array_t<double, py::array::f_style | py::array::forcecast> Q_arr,
    double min_mac,
    double min_var,
    int n_threads
) {
    return process_block_dense_impl(G_arr, rif_matrices, Q_arr, min_mac, min_var, n_threads);
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

    m.def("stage1_block_scores_multi", &stage1_block_scores_multi,
          "Compute block scores for multiple phenotypes",
          py::arg("G"), py::arg("rif_matrices"), py::arg("Q"),
          py::arg("block_ids"), py::arg("n_blocks"),
          py::arg("min_mac") = 5.0, py::arg("min_var") = 1e-8,
          py::arg("n_threads") = 1);

    m.def("process_block_dense_impl", &process_block_dense_impl,
          "Compute per-SNP betas/SEs for multiple phenotypes",
          py::arg("G"), py::arg("rif_matrices"), py::arg("Q"),
          py::arg("min_mac") = 5.0, py::arg("min_var") = 1e-8,
          py::arg("n_threads") = 1);

    m.def("process_block_sparse_impl", &process_block_sparse_impl,
          "Compute per-SNP betas/SEs for multiple phenotypes",
          py::arg("G"), py::arg("rif_matrices"), py::arg("Q"),
          py::arg("min_mac") = 5.0, py::arg("min_var") = 1e-8,
          py::arg("n_threads") = 1);
    
    m.def("openmp_info", &openmp_info, "Get OpenMP configuration");
}
