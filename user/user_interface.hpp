/* Copyright [2018] */

#ifndef INCLUDE_USER_INTERFACE_HPP_
#define INCLUDE_USER_INTERFACE_HPP_

#include <stdlib.h>
namespace user {

constexpr size_t opt_n = 4;            // number of optimization variables
constexpr size_t opt_m = 3;            // dim(Im(g))

template<class T>
void eval_obj(const T *x, T& obj_value) {
    obj_value = x[0] * x[0] + 2.0 * x[1] * x[1] - x[2];
}

template<class T>
void eval_constraints(const T *x, T* g) {
    g[0] = x[0] * x[0] + x[2] * x[2] + x[0] * x[2];
    g[1] = x[2] - x[3];
    g[2] = x[1] + x[3];
}

void get_starting_point(double* x, double* lam) {
    // initialize passive variables (ADOL-C)
    x[0]   = 2.0;
    x[1]   = 2.0;
    x[2]   = 1.0;
    x[3]   = 0.0;

    lam[0] = 0.0;
    lam[1] = 0.0;
    lam[2] = 0.0;
}

}  // namespace user

#endif  // INCLUDE_USER_INTERFACE_HPP_
