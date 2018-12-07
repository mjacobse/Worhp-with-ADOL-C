/* Copyright [2018] */

#ifndef CORE_GENERATE_TAPES_HPP_
#define CORE_GENERATE_TAPES_HPP_

#include <adolc/adolc.h>
#include <adolc/adolc_sparse.h>
#include <worhp/worhp.h>

#include "adolc_symbols.hpp"

// tamplated version of UserF that can process active and passive variables
template<class T>
void eval_obj(const T *x, T& obj_value);

// tamplated version of UserG that can process active and passive variables
template<class T>
void eval_constraints(const T *x, T *g);

void get_starting_point(double* x, double* lam);

void generate_tapes(int n, int m, int& nnz_grad_f, int& nnz_jac_g, int& nnz_h_lag, Workspace* wsp) {
    double *xp    = new double[n];
    double *lamp  = new double[m];
    double *zl    = new double[m];
    double *zu    = new double[m];

    adouble *xa   = new adouble[n];
    adouble *g    = new adouble[m];
    double *lam   = new double[m];
    double sig;
    adouble obj_value;

    double dummy;

    user::get_starting_point(xp, lamp);

    // taping the evaluation of the active counterpart to UserF
    trace_on(adolc::tag_f);
        // declare xa[i] as independent variables
        for (int idx = 0; idx < n; idx++)
            xa[idx] <<= xp[idx];

        user::eval_obj(xa, obj_value);

        // declare obj_value as dependent variable
        obj_value >>= dummy;
    trace_off();

    // taping the evaluation of the active counterpart to UserG
    trace_on(adolc::tag_g);
        for (int idx = 0; idx < n; idx++)
            xa[idx] <<= xp[idx];

        user::eval_constraints(xa, g);

        for (int idx = 0; idx < m; idx++)
            g[idx] >>= dummy;
    trace_off();

    // taping the evaluation of the active version of the Lagrangian
    trace_on(adolc::tag_L);
        for (int idx = 0; idx < n; idx++)
            xa[idx] <<= xp[idx];

        for (int idx = 0; idx < m; idx++)
            lam[idx] = lamp[idx];

        sig = wsp->ScaleObj;

        user::eval_obj(xa, obj_value);

        // explicit passive decalration of sig
        obj_value *= mkparam(sig);
        user::eval_constraints(xa, g);

        for (int idx = 0; idx < m; idx++)
            obj_value += g[idx]*mkparam(lam[idx]);

        obj_value >>= dummy;
    trace_off();

    adolc::rind_f = NULL;
    adolc::cind_f = NULL;

    adolc::rind_g = NULL;
    adolc::cind_g = NULL;

    adolc::rind_L = NULL;
    adolc::cind_L = NULL;

    adolc::gradval = NULL;
    adolc::jacval  = NULL;
    adolc::hessval = NULL;

    // computation of the sparsity pattern of gradient(UserF)
    sparse_jac(adolc::tag_f, 1, n, adolc::compute_pattern, xp, &adolc::nnz_grad_f,
              &adolc::rind_f, &adolc::cind_f, &adolc::gradval, adolc::options_g);
    wsp->DF.nnz = adolc::nnz_grad_f;

    // computation of the sparsity pattern of Jacobian(UserG)
    sparse_jac(adolc::tag_g, m, n, adolc::compute_pattern, xp, &adolc::nnz_jac_g,
              &adolc::rind_g, &adolc::cind_g, &adolc::jacval, adolc::options_g);
    wsp->DG.nnz = adolc::nnz_jac_g;

    // computation of the sparsity pattern of Hessian(Lagangian)
    sparse_hess(adolc::tag_L, n, adolc::compute_pattern, xp, &adolc::nnz_L,
               &adolc::rind_L, &adolc::cind_L, &adolc::hessval, adolc::options_L);

    // determine the additional sparsity etries of HM
    int additionalEntries4WorhpDiagonal = n;
    for (int i = 0; i < adolc::nnz_L; ++i)
        if (adolc::rind_L[i] == adolc::cind_L[i]) --additionalEntries4WorhpDiagonal;

    // output number of non-zeros in the hessian of the lagrangian
    adolc::nnz_h_lag = adolc::nnz_L + additionalEntries4WorhpDiagonal;
    wsp->HM.nnz      = adolc::nnz_h_lag;

    delete[] lam;
    delete[] g;
    delete[] xa;
    delete[] zu;
    delete[] zl;
    delete[] lamp;
    delete[] xp;
}

#endif  // CORE_GENERATE_TAPES_HPP_
