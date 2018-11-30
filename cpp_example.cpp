/*-----------------------------------------------------------------------
 *
 * Minimise    f
 *
 * subject to      -0.5 <= x1 <=  INFTY
 *                   -2 <= x2 <=  INFTY
 *                    0 <= x3 <=  2
 *                   -2 <= x4 <=  2
 *                         g1 ==  1
 *               -INFTY <= g2 <= -1
 *                  2.5 <= g3 <=  5
 *
 * where         f (x1,x2,x3,x4) = x1^2 + 2 x2^2 - x3
 *               g1(x1,x2,x3,x4) = x1^2 + x3^2 + x1x3
 *               g2(x1,x2,x3,x4) = x3 - x4
 *               g3(x1,x2,x3,x4) = x2 + x4
 *
 * Optimal solution
 *                     x*  = (0, 0.5, 1, 2)
 *                   f(x*) = -0.5
 *                   g(x*) = (1, -1, 2.5)
 *
 *-----------------------------------------------------------------------*/

#include "include/cpp_example.hpp" 
#include "include/sort_worhp_matrices.hpp"

int main()
{
    /*
	 * WORHP data structures
     *
     * OptVar contains the variables, constraint values and multipliers.
     * Workspace encapsulates all workspace and working variables for WORHP.
     * Params contains all WORHP parameters.
     * Control contains things for reverse communication flow control.
     */
    OptVar    opt;
    Workspace wsp;
    Params    par;
    Control   cnt;
 
    // Check Version of library and header files
    CHECK_WORHP_VERSION
 
    // Properly zeros everything, or else the following routines could get confused
    WorhpPreInit(&opt, &wsp, &par, &cnt);
 
    // Uncomment this to get more info on data structures
    //WorhpDiag(&opt, &wsp, &par, &cnt);
 
    /*
     * Parameter initialisation routine that must be called
     * when using ReadParamsNoInit instead of ReadParams.
     */
    int status;
    InitParams(&status, &par);
 
    /*
     * We can now set parameters that may be overruled by those in the
     * parameter file. This is useful for setting a non-default standard
     * parameter value that may still be overwritten.
     */
    par.NLPprint = 1;  // Let's prefer the slim output format
                       // unless the parameter file says differently
 
    /*
     * Parameter XML import routine that does not reset
     * all parameters to default values (InitParams does this)
     */
    ReadParamsNoInit(&status, "worhp.xml", &par);
    if (status == DataError || status == InitError)
    {
        return EXIT_FAILURE;
    }
 
    /*
     * WORHP data structure initialisation routine.
     * Calling this routine prior to WORHP is mandatory.
     * Before calling WorhpInit, set the problem and matrix dimensions as
     *
     * opt.n      = number of variables,
     * opt.m      = number of constraints (lin + nonlin, excluding box con's),
     * wsp.DF.nnz = nonzero entries of the objective function gradient,
     * wsp.DG.nnz = nonzero entries of the constraint Jacobian,
     * wsp.HM.nnz = nonzero entries of the Lagrange Hessian.
     *
     * Set nnz to 'WorhpMatrix_Init_Dense' to have WorhpInit allocate and
     * create a dense matrix structure appropriate for the matrix kind and
     * its dimensions. Setting it to its dense dimension achieves the same.
     */
    opt.n = 4;  // This problem has 4 variables
    opt.m = 3;  // and 3 constraints (excluding box constraints)
 
    /*
     * ADOLC data structure initialisation raoutine.
     * Calling this routine prior to any exploitation of automatic derivatives is mandatory.
     * Input to generate_tapes: are the dimensions of the problem (dim(UserF^-1) and dim(Im(UserG)))
     * I/O containers storing the number of non-zeros in the jacobian and the hessian 
     * and a poiter to the WORHP workspace.
     * 
     * NOTE: nnz_h_lag is already considers the full diagonal indepented of its sparsity pattern
     */ 

    generate_tapes(opt.n, opt.m, nnz_jac_g, nnz_h_lag, &wsp);
 
    wsp.DF.nnz = 3;
    wsp.DG.nnz = nnz_jac_g;
    wsp.HM.nnz = nnz_h_lag;  

    WorhpInit(&opt, &wsp, &par, &cnt);
    if (cnt.status != FirstCall)
    {
        std::cout << "Main: Initialisation failed." << std::endl;
        return EXIT_FAILURE;
    }
 
    /*
     * These pointers give access to the essential user data:
     *
     * opt.X[0] to opt.X[opt.n - 1]           : Optimisation variables
     * opt.Lambda[0] to opt.Lambda[opt.n - 1] : Multipliers for the constraints
     *                                          on X ("box constraints")
     * opt.G[0] to opt.G[opt.m - 1]           : Linear and nonlinear constraints
     * opt.Mu[0] to opt.Mu[opt.m - 1]         : Multipliers for the constraints on G
     *
     * Set initial values of X, Lambda and Mu here.
     * G need not be initialised.
     */
    opt.X[0] = 2.0;
    opt.X[1] = 2.0;
    opt.X[2] = 1.0;
    opt.X[3] = 0.0;
    opt.Lambda[0] = 0.0;
    opt.Lambda[1] = 0.0;
    opt.Lambda[2] = 0.0;
    opt.Lambda[3] = 0.0;
    opt.Mu[0] = 0.0;
    opt.Mu[1] = 0.0;
    opt.Mu[2] = 0.0;
 
    /*
     * Set lower and upper bounds on the variables and constraints.
     * Use +/-par.Infty to signal "unbounded".
     *
     * XL and XU are lower and upper bounds ("box constraints") on X.
     * GL and GU are lower and upper bounds on G.
     */
    opt.XL[0] = -0.5;
    opt.XU[0] =  par.Infty;
    opt.XL[1] = -2.0;
    opt.XU[1] =  par.Infty;
    opt.XL[2] =  0.0;
    opt.XU[2] =  2.0;
    opt.XL[3] = -2.0;
    opt.XU[3] =  2.0;
 
    opt.GL[0] =  1.0;  // set opt.GL[i] == opt.GU[i]
    opt.GU[0] =  1.0;  // for equality constraints
    opt.GL[1] = -par.Infty;
    opt.GU[1] = -1.0;
    opt.GL[2] =  2.5;
    opt.GU[2] =  5.0;
 
    /*
     * Specify matrix structures in CS format, using Fortran indexing,
     * i.e. 1...N instead of 0...N-1, to describe the matrix structure.
     */
 
    // Define DF as row vector, column index is ommited
    if (wsp.DF.NeedStructure)
    {
        // only set the nonzero entries, so omit the 4th entry,
        // which is a structural zero
        wsp.DF.row[0] = 1;
        wsp.DF.row[1] = 2;
        wsp.DF.row[2] = 3;
    }
 
    // Define DG as CS-matrix
    if (wsp.DG.NeedStructure)
    {
        // only set the nonzero entries in column-major order
        
        /*
         * This piece of code is intended to reorder the arrays rind_g and cind_g. ADOL-C
         * returns them in row major order WORHP needs them in column order. The three components
         * of the sparsity pattern rowIdx (rind_g), colIdx (cind_g) and value at this position (jacval).
         * They are gathered in tupels these tupel are then ordered in column-major order.
         * Finally the sparsity pattern of DG is filled with with values.
         */ 
        
        std::vector<MatrixEntry> sparseDG;
        for(int i =0; i <nnz_jac_g; ++i)
        {
            sparseDG.emplace_back(Row(rind_g[i]), Col(cind_g[i]), Value(jacval[i]));
        }

        std::sort(sparseDG.begin(),sparseDG.end(), sortDG<>());

        for(int i =0; i <nnz_jac_g; ++i)
        {
            wsp.DG.row[i] = sparseDG[i].getRow().to_int() +1;
            wsp.DG.col[i] = sparseDG[i].getCol().to_int() +1;
        }
    }
 
    // Define HM as a diagonal LT-CS-matrix, but only if needed by WORHP
    if (wsp.HM.NeedStructure)
    {
        /*
         * This piece of code is intended to reorder the arrays rind_L and cind_L. ADOL-C
         * returns them in row-major order and the upper triangular matrix.
         * WORHP needs them in a special order. Each diagonal element is bigger than any
         * non-diagonal element. Internally both partitions are sorted in column-major order.
         * 
         * NOTE: WORHP requires a full diagonal independent of this actual sparsity structure,
         * this is considered in this block of code. After the elements are sorted they written to
         * The wsp.HM pattern with tranposed indices.
         */ 

        std::vector<MatrixEntry> sparseHM;
        std::set<int> missingDiagonalElems {};

        for(int i = 0; i < opt.n; ++i)
        {
            missingDiagonalElems.insert(i);
        }

        for(int i =0; i < nnz_L; ++i)
        {
            sparseHM.emplace_back(Row(rind_L[i]), Col(cind_L[i]), Value(hessval[i]));
            
            if(rind_L[i] == cind_L[i])
            {
                missingDiagonalElems.erase(rind_L[i]);
            }
        }

        for(const auto idx : missingDiagonalElems)
        {
            sparseHM.emplace_back(Row(idx), Col(idx), Value(0));
        }

        std::sort(sparseHM.begin(),sparseHM.end(), sortHM<>(opt.n));

        for(size_t i = 0; i < sparseHM.size(); ++i)
        {
            wsp.HM.row[i] = sparseHM[i].getCol().to_int() + 1;
            wsp.HM.col[i] = sparseHM[i].getRow().to_int() + 1;
        }
    }
 
    /*
     * WORHP Reverse Communication loop.
     * In every iteration poll GetUserAction for the requested action, i.e. one
     * of {callWorhp, iterOutput, evalF, evalG, evalDF, evalDG, evalHM, fidif}.
     *
     * Make sure to reset the requested user action afterwards by calling
     * DoneUserAction, except for 'callWorhp' and 'fidif'.
     */
    while (cnt.status < TerminateSuccess && cnt.status > TerminateError)
    {
        /*
         * WORHP's main routine.
         * Do not manually reset callWorhp, this is only done by the FD routines.
         */
        if (GetUserAction(&cnt, callWorhp))
        {
            Worhp(&opt, &wsp, &par, &cnt);
            // No DoneUserAction!
        }
 
        /*
         * Show iteration output.
         * The call to IterationOutput() may be replaced by user-defined code.
         */
        if (GetUserAction(&cnt, iterOutput))
        {
            IterationOutput(&opt, &wsp, &par, &cnt);
            DoneUserAction(&cnt, iterOutput);
        }
 
        /*
         * Evaluate the objective function.
         * The call to UserF may be replaced by user-defined code.
         */
        if (GetUserAction(&cnt, evalF))
        {
            UserF(&opt, &wsp, &par, &cnt);
            DoneUserAction(&cnt, evalF);
        }
 
        /*
         * Evaluate the constraints.
         * The call to UserG may be replaced by user-defined code.
         */
        if (GetUserAction(&cnt, evalG))
        {
            UserG(&opt, &wsp, &par, &cnt);
            DoneUserAction(&cnt, evalG);
        }
 
        /*
         * Evaluate the gradient of the objective function.
         * The call to UserDF may be replaced by user-defined code.
         */
        if (GetUserAction(&cnt, evalDF))
        {
            UserDF(&opt, &wsp, &par, &cnt);
            DoneUserAction(&cnt, evalDF);
        }
 
        /*
         * Evaluate the Jacobian of the constraints.
         * The call to UserDG may be replaced by user-defined code.
         */
        if (GetUserAction(&cnt, evalDG))
        {
            UserDG(&opt, &wsp, &par, &cnt);
            DoneUserAction(&cnt, evalDG);
        }
 
        /*
         * Evaluate the Hessian matrix of the Lagrange function (L = f + mu*g)
         * The call to UserHM may be replaced by user-defined code.
         */
        if (GetUserAction(&cnt, evalHM))
        {
            UserHM(&opt, &wsp, &par, &cnt);
            DoneUserAction(&cnt, evalHM);
        }
 
        /*
         * Use finite differences with RC to determine derivatives
         * Do not reset fidif, this is done by the FD routine.
         */
        if (GetUserAction(&cnt, fidif))
        {
            WorhpFidif(&opt, &wsp, &par, &cnt);
            // No DoneUserAction!
        }
    }
 
    // Translate the WORHP status flag into a meaningful message.
    StatusMsg(&opt, &wsp, &par, &cnt);
    // Deallocate all data structures.
    // Data structures must not be accessed after this call.
    WorhpFree(&opt, &wsp, &par, &cnt);
 
    return EXIT_SUCCESS;
}
 

template<class T>
bool eval_obj(const T *x, T& obj_value)
{
	obj_value = (x[0] * x[0] + 2.0 * x[1] * x[1] - x[2]);
    return true;
}

template<class T>
bool eval_constraints(const T *x, T* g)
{
    g[0] = x[0] * x[0] + x[2] * x[2] + x[0] * x[2];
    g[1] = x[2] - x[3];
    g[2] = x[1] + x[3];
    return true;
}

void UserF(OptVar *opt, Workspace *wsp, Params *par, Control *cnt)
{
    double *X = opt->X;
    double obj_value;
    eval_obj(X,obj_value);
    opt->F = wsp->ScaleObj * obj_value;
}
 
void UserG(OptVar *opt, Workspace *wsp, Params *par, Control *cnt)
{
    double *X = opt->X;
    double g[3];
    
    eval_constraints(X,g);

    opt->G[0] = g[0];
    opt->G[1] = g[1];
    opt->G[2] = g[2];
}
 
void UserDF(OptVar *opt, Workspace *wsp, Params *par, Control *cnt)
{
    wsp->DF.val[0] = wsp->ScaleObj *  2.0 * opt->X[0];
    wsp->DF.val[1] = wsp->ScaleObj *  4.0 * opt->X[1];
    wsp->DF.val[2] = wsp->ScaleObj * -1.0;
}
 
void UserDG(OptVar *opt, Workspace *wsp, Params *par, Control *cnt)
{
    double *X = opt->X;

    sparse_jac(tag_g, opt->m, opt->n, reuse_pattern, X, &nnz_jac, &rind_g, &cind_g, &jacval, options_g); 

    
    std::vector<MatrixEntry> sparseDG;
    for(int i =0; i <nnz_jac_g; ++i)
    {
        sparseDG.emplace_back(Row(rind_g[i]), Col(cind_g[i]), Value(jacval[i]));    
    }

    std::sort(sparseDG.begin(),sparseDG.end(), sortDG<>());
        
    for(int i =0; i <nnz_jac_g; ++i)
    {
        wsp->DG.val[i] = sparseDG[i].getVal().to_double();
    }
}
 
void UserHM(OptVar *opt, Workspace *wsp, Params *par, Control *cnt)
{
    double *X = opt->X;  

    sparse_hess(tag_L, opt->n, reuse_pattern, X, &nnz_L, &rind_L, &cind_L, &hessval, options_L);

    std::vector<MatrixEntry>  sparseHM;

    std::set<int> missingDiagonalElems {};
    for(int i = 0; i < opt->n; ++i )
    {
        missingDiagonalElems.insert(i);
    }

    for(int i =0; i < nnz_L; ++i)
    {
        sparseHM.emplace_back(Row(rind_L[i]), Col(cind_L[i]), Value(hessval[i]));
        
        if(rind_L[i] == cind_L[i])
        {
            missingDiagonalElems.erase(rind_L[i]);
        }
    }

    for(const auto idx : missingDiagonalElems)
    {
        sparseHM.emplace_back(Row(idx),Col(idx),Value(0));
    }

    std::sort(sparseHM.begin(),sparseHM.end(), sortHM<>(opt->n));

    for(size_t i =0; i < sparseHM.size(); ++i)
    {
            wsp->HM.val[i] = sparseHM[i].getVal().to_double();
    }
}


void generate_tapes(int n, int m, int& nnz_jac_g, int& nnz_h_lag, Workspace* wsp)
{
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

    // initialize passive variables
    xp[0]   = 2.0;
    xp[1]   = 2.0;
    xp[2]   = 1.0;
    xp[3]   = 0.0;

    lamp[0] = 0.0;
    lamp[1] = 0.0;
    lamp[2] = 0.0;


    // taping the evaluation of the active counterpart to UserF
    trace_on(tag_f);
        // declare xa[i] as independent variables
        for(int idx=0; idx<n; idx++)
            xa[idx] <<= xp[idx];

        eval_obj(xa,obj_value);
        
        // declare obj_value as dependent variable
        obj_value >>= dummy;
    trace_off();
    
    // taping the evaluation of the active counterpart to UserG
    trace_on(tag_g);
        for(int idx=0;idx<n;idx++)
            xa[idx] <<= xp[idx];

        eval_constraints(xa,g);

        for(int idx=0;idx<m;idx++)
            g[idx] >>= dummy;
    trace_off();

    // taping the evaluation of the active version of the Lagrangian
    trace_on(tag_L); 
       for(int idx=0;idx<n;idx++)
            xa[idx] <<= xp[idx];
    
        for(int idx=0;idx<m;idx++)
            lam[idx] = lamp[idx];
    
        sig = wsp->ScaleObj;

        eval_obj(xa,obj_value);
        
        // explicit passive decalration of sig
        obj_value *= mkparam(sig);
        eval_constraints(xa,g);
 
        for(int idx=0;idx<m;idx++)
            obj_value += g[idx]*mkparam(lam[idx]);

        obj_value >>= dummy;
    trace_off();

    rind_g = NULL; 
    cind_g = NULL;
    rind_L = NULL;
    cind_L = NULL;

    jacval  = NULL;
    hessval = NULL;
  
    // computation of the sparsity pattern of Jacobian(UserG)
    sparse_jac(tag_g, m, n, compute_pattern, xp, &nnz_jac, &rind_g, &cind_g, &jacval, options_g); 

    // output number of non-zeros in the jacobian of userG
    nnz_jac_g = nnz_jac;

    // computation of the sparsity pattern of Hessian(Lagangian)  
    sparse_hess(tag_L, n, compute_pattern, xp, &nnz_L, &rind_L, &cind_L, &hessval, options_L);


    // determine the additional sparsity etries of HM
    int additionalEntries4Worhp = n;
    for(int i=0; i < nnz_L; ++i)
    {
        if(rind_L[i] == cind_L[i]) --additionalEntries4Worhp;
    }
    // output number of non-zeros in the hessian of the lagrangian
    nnz_h_lag = nnz_L + additionalEntries4Worhp;

    delete[] lam;
    delete[] g;
    delete[] xa;
    delete[] zu;
    delete[] zl;
    delete[] lamp;
    delete[] xp;
}