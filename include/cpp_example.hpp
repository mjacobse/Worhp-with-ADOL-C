#ifndef CPP_EXAMPLE_HPP
  #define CPP_EXAMPLE_HPP
  
  #include <iostream>
  #include <array>
  #include <algorithm>
  #include <tuple>
  #include <unordered_set>
  #include <set>
  #include <cassert>
  #include <cstring>
  #include "worhp/worhp.h"
  #include <adolc/adolc.h>
  #include <adolc/adolc_sparse.h>

  // defines for adolc constants
  constexpr int tag_f = 1;            // identifier for the tape that stores the objective function UserF
  constexpr int tag_g = 2;            // identifier for the tape that stores the constraints function UserG
  constexpr int tag_L = 3;            // identifier for the tape that stores the Lagrangian

  constexpr int compute_pattern = 0;  // indicates that the sparsity pattern shall be computed
  constexpr int reuse_pattern   = 1;  // reuse that the sparsity pattern shall be computed


  //** variables for sparsity exploitation
  unsigned int *rind_g;        // Jacobian g: row indices jacobian
  unsigned int *cind_g;        // Jacobian g: column indices
  double *jacval;              // Jacobian g: values

  unsigned int *rind_L;        // Hessian L: row indices
  unsigned int *cind_L;        // Hessian L: column indices
  double *hessval;             // Hessian L: values

  int nnz_jac_g;               // number of non zeros in the jacobian of UserG
  int nnz_h_lag;               // number of non zeros in lower triangular of the hessian of the lagrangian

  int nnz_jac;                 // number of non zeros in the jacobian of UserG (used only in alocal scope, can be removed eventually )
  int nnz_L;                   // number of non zeros in lower triangular of the hessian of the lagrangian (used only in alocal scope, can be removed eventually )


  int options_g[4] = {0,0,0,0}; // parameter vector for sparse_jac
  int options_L[2] = {0,1};     // parameter vector for sparse_hessian


  // tamplated version of UserF that can process active and passive variables
  template<class T>
  bool eval_obj(const T *x, T& obj_value);

  // tamplated version of UserG that can process active and passive variables
  template<class T>
  bool eval_constraints(const T *x, T *g);

  // genretes the tapes for the functions UserF, UserG and the corresponding Lagrangian
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

#endif