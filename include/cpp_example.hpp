/* Copyright (C) 2018 Luis Lüttgens - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the EPL or GPL license.
 *
 * You should have received a copy of the XYZ license with
 * this file. If not, please write to: luis.luett@googlemail.com
 */

#ifndef INCLUDE_CPP_EXAMPLE_HPP_
    #define INCLUDE_CPP_EXAMPLE_HPP_

    // .h header
    #include <adolc/adolc.h>
    #include <adolc/adolc_sparse.h>
    #include <worhp/worhp.h>

    // cpp header

    #include <cassert>
    #include <cstring>

    #include <iostream>
    #include <array>
    #include <algorithm>
    #include <tuple>
    #include <unordered_set>
    #include <set>



    // defines for adolc constants
    constexpr int tag_f = 1;            // tape identifier for UserF
    constexpr int tag_g = 2;            // tape identifier for UserG
    constexpr int tag_L = 3;            // tape identifier for UserHM

    constexpr int compute_pattern = 0;  // computed the sparsity pattern
    constexpr int reuse_pattern   = 1;  // reuse the sparsity pattern


    //** variables for sparsity exploitation
    unsigned int *rind_g;        // Jacobian g: row indices jacobian
    unsigned int *cind_g;        // Jacobian g: column indices
    double *jacval;              // Jacobian g: values

    unsigned int *rind_L;        // Hessian L: row indices
    unsigned int *cind_L;        // Hessian L: column indices
    double *hessval;             // Hessian L: values

    int nnz_jac_g;               // number of non zeros in UserDG
    int nnz_h_lag;               // number of non zeros in UserHM

    int nnz_jac;                 // number of non zeros in UserDG (local scope)
    int nnz_L;                   // number of non zeros in UserHM (local scope)


    int options_g[4] = {0, 0, 0, 0};  // parameter vector for sparse_jac
    int options_L[2] = {0, 1};        // parameter vector for sparse_hessian


    // tamplated version of UserF that can process active and passive variables
    template<class T>
    bool eval_obj(const T *x, T& obj_value);

    // tamplated version of UserG that can process active and passive variables
    template<class T>
    bool eval_constraints(const T *x, T *g);

    // genretes the tapes tag_f, tag_g and tag_L
    void generate_tapes(int n, int m, int &nnz_jac_g, int &nnz_h_lag, Workspace* wsp);

    // Objective function
    void UserF(OptVar *opt, Workspace *wsp, Params *par, Control *cnt);

    // Function of constraints
    void UserG(OptVar *opt, Workspace *wsp, Params *par, Control *cnt);

    // Gradient of objective function
    void UserDF(OptVar *opt, Workspace *wsp, Params *par, Control *cnt);

    // Jacobian of constraints
    void UserDG(OptVar *opt, Workspace *wsp, Params *par, Control *cnt);

    // Hessian of Lagrangian
    void UserHM(OptVar *opt, Workspace *wsp, Params *par, Control *cnt);

#endif  // INCLUDE_CPP_EXAMPLE_HPP_
