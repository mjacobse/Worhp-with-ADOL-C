/* Copyright [2018] */
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

#include "cpp_example.hpp"
#include "include/worhpAD.hpp"

int main() {
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

    // Properly zeros everything, or else the following
    //  routines could get confused

    WorhpPreInit(&opt, &wsp, &par, &cnt);

    // Uncomment this to get more info on data structures
    // WorhpDiag(&opt, &wsp, &par, &cnt);

    /*
     * Parameter initialisation routine that must be called
     * when using ReadParamsNoInit instead of ReadParams.
     */

    int status {};
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
    if (status == DataError || status == InitError) {
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
    opt.n = user::opt_n;  // This problem has 4 variables
    opt.m = user::opt_m;  // and 3 constraints (excluding box constraints)

    /*
     * ADOLC data structure initialisation raoutine.
     * Calling this routine prior to any exploitation of automatic derivatives is mandatory.
     * Input to generate_tapes: are the dimensions of the problem (dim(UserF^-1) and dim(Im(UserG)))
     * I/O containers storing the number of non-zeros in the jacobian and the hessian 
     * and a poiter to the WORHP workspace.
     * 
     * NOTE: nnz_h_lag is already considers the full diagonal indepented of its sparsity pattern
     */ 

    generate_tapes(user::opt_n, user::opt_m, adolc::nnz_grad_f,
                   adolc::nnz_jac_g, adolc::nnz_h_lag, &wsp);

    WorhpInit(&opt, &wsp, &par, &cnt);
    if (cnt.status != FirstCall) {
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
    if (wsp.DF.NeedStructure) {
        worhp::auto_diff_DF_pattern(&wsp);
    }

    // Define DG as CS-matrix
    if (wsp.DG.NeedStructure) {
        worhp::auto_diff_DG_pattern(&wsp);
    }

    // Define HM as a diagonal LT-CS-matrix, but only if needed by WORHP
    if (wsp.HM.NeedStructure) {
        worhp::auto_diff_HM_pattern(&wsp);
    }

    /*
     * WORHP Reverse Communication loop.
     * In every iteration poll GetUserAction for the requested action, i.e. one
     * of {callWorhp, iterOutput, evalF, evalG, evalDF, evalDG, evalHM, fidif}.
     *
     * Make sure to reset the requested user action afterwards by calling
     * DoneUserAction, except for 'callWorhp' and 'fidif'.
     */
    while (cnt.status < TerminateSuccess && cnt.status > TerminateError) {
        /*
         * WORHP's main routine.
         * Do not manually reset callWorhp, this is only done by the FD routines.
         */
        if (GetUserAction(&cnt, callWorhp)) {
            Worhp(&opt, &wsp, &par, &cnt);
            // No DoneUserAction!
        }

        /*
         * Show iteration output.
         * The call to IterationOutput() may be replaced by user-defined code.
         */
        if (GetUserAction(&cnt, iterOutput)) {
            IterationOutput(&opt, &wsp, &par, &cnt);
            DoneUserAction(&cnt, iterOutput);
        }

        /*
         * Evaluate the objective function.
         * The call to UserF may be replaced by user-defined code.
         */
        if (GetUserAction(&cnt, evalF)) {
            UserF(&opt, &wsp, &par, &cnt);
            DoneUserAction(&cnt, evalF);
        }

        /*
         * Evaluate the constraints.
         * The call to UserG may be replaced by user-defined code.
         */
        if (GetUserAction(&cnt, evalG)) {
            UserG(&opt, &wsp, &par, &cnt);
            DoneUserAction(&cnt, evalG);
        }

        /*
         * Evaluate the gradient of the objective function.
         * The call to UserDF may be replaced by user-defined code.
         */
        if (GetUserAction(&cnt, evalDF)) {
            UserDF(&opt, &wsp, &par, &cnt);
            DoneUserAction(&cnt, evalDF);
        }

        /*
         * Evaluate the Jacobian of the constraints.
         * The call to UserDG may be replaced by user-defined code.
         */
        if (GetUserAction(&cnt, evalDG)) {
            UserDG(&opt, &wsp, &par, &cnt);
            DoneUserAction(&cnt, evalDG);
        }

        /*
         * Evaluate the Hessian matrix of the Lagrange function (L = f + mu*g)
         * The call to UserHM may be replaced by user-defined code.
         */
        if (GetUserAction(&cnt, evalHM)) {
            UserHM(&opt, &wsp, &par, &cnt);
            DoneUserAction(&cnt, evalHM);
        }

        /*
         * Use finite differences with RC to determine derivatives
         * Do not reset fidif, this is done by the FD routine.
         */
        if (GetUserAction(&cnt, fidif)) {
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
