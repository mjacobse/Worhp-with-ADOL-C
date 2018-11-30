This project presents a way to combine the NLP solver WORHP and the automatic differentiation tool ADOL-C. We have reimplemented the "Getting started" example of WOHRP available [here](https://worhp.de/content/cppexample), but instead of computing the required derivatives by hand, we are using ADOL-C drivers. ADOL-C is available [here](https://www.coin-or.org/download/source/ADOL-C/)

# Motivation
The used case of this combination is rapid prototyping. Modeling a dynamic process is usually an iterative process. The complexity of the model rises gradually. We could either recompute all derivatives analytically every time the model changes or use ADOL-C to do this automatically for us.

# Code Style
This project tries to follow the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html) enforced by  [cpplint](https://github.com/cpplint/cpplint).

# Code Example
ADOL-C makes use of active (`adouble`) and passive data types (any other type). Integrating ADOL-C into  WORHP code (or any other NLP solver for that matter) consists of three steps.

## Overload UserF and UserG
ADOL-C monitors all differentiable functions on independent tapes. For that purpose, all of the functions are evaluated with the adoubles. During the optimization, the same functions are being evaluated with their passive counterparts. For that reason, it makes sense to use template functions. But we still have to leave the interface to WORHP intact. We outsource the computations to an independent function and calling inside WORHP functions these functions. This has the following effect one `UserF`.
```
void UserF(OptVar *opt, Workspace *wsp, Params *par, Control *cnt)
{
    double *X = opt->X;  // Abbreviate notation
 
    opt->F = wsp->ScaleObj * (X[0] * X[0] + 2.0 * X[1] * X[1] - X[2]);
}
```
Becomes the following block of code.
```
template<class T>
bool eval_obj(const T *x, T& obj_value) {
    obj_value = (x[0] * x[0] + 2.0 * x[1] * x[1] - x[2]);
    return true;
}

void UserF(OptVar *opt, Workspace *wsp, Params *par, Control *cnt) {
    double *X = opt->X;
    double obj_value;
    eval_obj(X, obj_value);
    opt->F = wsp->ScaleObj * obj_value;
}
```

## Generation of the tapes
In this step we are generating the tapes for our three functions of interest:  `UserF, UserG` and Lagrangian. We activate the monitoring and then evaluate each of this functions with active variables.
```
    // taping the evaluation of the active counterpart to UserG
    trace_on(tag_g);
        for (int idx = 0; idx < n; idx++)
            xa[idx] <<= xp[idx];

        eval_constraints(xa, g);

        for (int idx = 0; idx < m; idx++)
            g[idx] >>= dummy;
    trace_off()
```
This procedure is not limited to a single function evaluation. Any combination of functions can be monitored (see `tag_L`). We are using the monitored information right after. With the help of ADOL-C's (sparse) drivers for optimization, we can compute either the Jacobian or the Hessian of the taped function.
```
// computation of the sparsity pattern of Jacobian(UserG)
    sparse_jac(tag_g, m, n, compute_pattern, xp, &nnz_jac,
                      &rind_g, &cind_g, &jacval, options_g);
```
ADOL-C stores the sparsity pattern in three independet arrays:
`rind_g, cind_g` and `jacval`

 ## Use Jacobian and Hessian
 
Now we want to use the Jacobian and Hessian inside WORHP. This is a combination of the previous two steps, instead of doing the evaluation by hand we are calling the same functions `sparse_jac` and `sparse_hess` as in the previous chapter. Their results are stored in an data structure called `MatrixEntry`. This data structure is used to sort the array `rind_L, cind_L` and `hessval` in the order WORHP requires. they are then written into at appropriate positions in `UserDG` and `UserHM`.
 
 # License
 
 This project is under the GPL license, see the LICENSE file for more information.
