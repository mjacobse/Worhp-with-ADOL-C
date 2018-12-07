/* Copyright [2018] */
#ifndef INCLUDE_ADOLC_SYMBOLS_HPP_
#define INCLUDE_ADOLC_SYMBOLS_HPP_
namespace adolc {
    // defines for adolc constants
    constexpr int tag_f = 1;            // tape identifier for UserF
    constexpr int tag_g = 2;            // tape identifier for UserG
    constexpr int tag_L = 3;            // tape identifier for UserHM

    constexpr int compute_pattern = 0;  // computed the sparsity pattern
    constexpr int reuse_pattern   = 1;  // reuse the sparsity pattern

    //** variables for sparsity exploitation
    unsigned int *rind_f;        // gradient f: row indices gradient
    unsigned int *cind_f;        // gradient f: column indices
    double *gradval;             // gradient f: values

    unsigned int *rind_g;        // Jacobian g: row indices jacobian
    unsigned int *cind_g;        // Jacobian g: column indices
    double *jacval;              // Jacobian g: values

    unsigned int *rind_L;        // Hessian L: row indices
    unsigned int *cind_L;        // Hessian L: column indices
    double *hessval;             // Hessian L: values

    int nnz_grad_f;              // number of non zeros in UserDF
    int nnz_jac_g;               // number of non zeros in UserDG
    int nnz_h_lag;               // augmented number of non zeros in UserHM (inclduing diagonal)

    int nnz_L;                   // actual number of non zeros in UserHM


    int options_g[4] = {0, 0, 0, 0};  // parameter vector for sparse_jac
    int options_L[2] = {0, 1};        // parameter vector for sparse_hessian

}  // namespace adolc

#endif  // INCLUDE_ADOLC_SYMBOLS_HPP_
