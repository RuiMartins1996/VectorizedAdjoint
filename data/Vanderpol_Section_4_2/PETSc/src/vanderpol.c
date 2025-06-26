static char help[] = "Performs adjoint sensitivity analysis for the van der Pol equation.\n";

/* ------------------------------------------------------------------------

   This program solves the van der Pol DAE ODE equivalent
      [ u_1' ] = [          u_2                ]  (2)
      [ u_2' ]   [ \mu ((1 - u_1^2) u_2 - u_1) ]
   on the domain 0 <= x <= 1, with the boundary conditions
       u_1(0) = 2, u_2(0) = - 2/3 +10/(81*\mu) - 292/(2187*\mu^2),
   and
       \mu = 10^6 ( y'(0) ~ -0.6666665432100101).,
   and computes the sensitivities of the final solution w.r.t. initial conditions and parameter \mu with the implicit theta method and its discrete adjoint.

   In an IMEX setting,  we can split (2) by component,

   [ u_1' ] = [ u_2 ] + [            0                ]
   [ u_2' ]   [  0  ]   [ \mu ((1 - u_1^2) u_2 - u_1) ]

   where

   [ G(u,t) ] = [ u_2 ]
                [  0  ]

   and

   [ F(u',u,t) ] = [ u_1' ] - [            0                ]
                   [ u_2' ]   [ \mu ((1 - u_1^2) u_2 - u_1) ]

   Notes:
   This code demonstrates the TSAdjoint interface to a DAE system.

   The user provides the implicit right-hand-side function
   [ F(u',u,t) ] = [u' - f(u,t)] = [ u_1'] - [        u_2             ]
                                   [ u_2']   [ \mu ((1-u_1^2)u_2-u_1) ]

   and the Jacobian of F (from the PETSc user manual)

              dF   dF
   J(F) = a * -- + --
              du'  du

   and the JacobianP of the explicit right-hand side of (2) f(u,t) ( which is equivalent to -F(0,u,t)).
   df   [       0               ]
   -- = [                       ]
   dp   [ (1 - u_1^2) u_2 - u_1 ].

   See ex20.c for more details on the Jacobian.

  ------------------------------------------------------------------------- */
#include <petscts.h>
#include <petsctao.h>



extern PetscErrorCode PetscLogView_VecScatter(PetscViewer);

typedef struct _n_User *User;
struct _n_User {
  PetscReal mu;
  PetscReal next_output;
  PetscBool imex;
  /* Sensitivity analysis support */
  PetscInt  steps;
  PetscReal ftime;
  Mat       A;                    /* IJacobian matrix */
  Mat       B;                    /* RHSJacobian matrix */
  Mat       Jacp;                 /* IJacobianP matrix */
  Mat       Jacprhs;              /* RHSJacobianP matrix */
  Vec       U, lambda[2], mup[2]; /* adjoint variables */
  PetscReal max_dt;
};

/* ----------------------- Explicit form of the ODE  -------------------- */

static PetscErrorCode RHSFunction(TS ts, PetscReal t, Vec U, Vec F, void *ctx)
{
  User               user = (User)ctx;
  PetscScalar       *f;
  const PetscScalar *u;

  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(U, &u));
  PetscCall(VecGetArray(F, &f));
  f[0] = u[1];
  if (user->imex) {
    f[1] = 0.0;
  } else {
    f[1] = ((1. - u[0] * u[0]) * u[1] - u[0])/user->mu;
  }
  PetscCall(VecRestoreArrayRead(U, &u));
  PetscCall(VecRestoreArray(F, &f));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RHSJacobian(TS ts, PetscReal t, Vec U, Mat A, Mat B, void *ctx)
{
  User               user     = (User)ctx;
  PetscReal          mu       = user->mu;
  PetscInt           rowcol[] = {0, 1};
  PetscScalar        J[2][2];
  const PetscScalar *u;

  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(U, &u));
  J[0][0] = 0;
  J[0][1] = 1.0;
  if (user->imex) {
    J[1][0] = 0.0;
    J[1][1] = 0.0;
  } else {
    J[1][0] = -(2.0 * u[1] * u[0] + 1.)/mu;
    J[1][1] =  (1.0 - u[0] * u[0])/mu;
  }
  PetscCall(MatSetValues(A, 2, rowcol, 2, rowcol, &J[0][0], INSERT_VALUES));
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  if (A != B) {
    PetscCall(MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY));
  }
  PetscCall(VecRestoreArrayRead(U, &u));
  PetscFunctionReturn(PETSC_SUCCESS);
}



/* ----------------------- Implicit form of the ODE  -------------------- */

static PetscErrorCode IFunction(TS ts, PetscReal t, Vec U, Vec Udot, Vec F, void *ctx)
{
  User               user = (User)ctx;
  const PetscScalar *u, *udot;
  PetscScalar       *f;
  PetscReal          mu       = user->mu;

  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(U, &u));
  PetscCall(VecGetArrayRead(Udot, &udot));
  PetscCall(VecGetArray(F, &f));
  if (user->imex) {
    f[0] = udot[0];
  } else {
    f[0] = udot[0] - u[1];
  }
  f[1] = udot[1] - ((1.0 - u[0] * u[0]) * u[1] - u[0])/mu;
  PetscCall(VecRestoreArrayRead(U, &u));
  PetscCall(VecRestoreArrayRead(Udot, &udot));
  PetscCall(VecRestoreArray(F, &f));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode IJacobian(TS ts, PetscReal t, Vec U, Vec Udot, PetscReal a, Mat A, Mat B, void *ctx)
{
  User               user     = (User)ctx;
  PetscInt           rowcol[] = {0, 1};
  PetscScalar        J[2][2];
  const PetscScalar *u;
  PetscReal          mu       = user->mu;


  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(U, &u));

  if (user->imex) {
    J[0][0] = a;
    J[0][1] = 0.0;
  } else {
    J[0][0] = a;
    J[0][1] = -1.0;
  }
  J[1][0] = (2.0 * u[0] * u[1] + 1.0)/mu;
  J[1][1] = a - (1.0 - u[0] * u[0])/mu;

  PetscCall(MatSetValues(B, 2, rowcol, 2, rowcol, &J[0][0], INSERT_VALUES));
  PetscCall(VecRestoreArrayRead(U, &u));

  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  if (B && A != B) {
    PetscCall(MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode IJacobianP(TS ts, PetscReal t, Vec U, Vec Udot, PetscReal a, Mat A, void *ctx)
{
  User               user     = (User)ctx;
  PetscInt           row[] = {0, 1}, col[] = {0};
  PetscScalar        J[2][1];
  const PetscScalar *u;

  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(U, &u));
  J[0][0] = 0;
  J[1][0] = ((1.0 - u[0] * u[0]) * u[1] - u[0])/(user->mu * user->mu);
  PetscCall(MatSetValues(A, 2, row, 1, col, &J[0][0], INSERT_VALUES));
  PetscCall(VecRestoreArrayRead(U, &u));
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RHSJacobianP(TS ts, PetscReal t, Vec U, Mat A, void *ctx)
{
  User user = (User)ctx;

  PetscFunctionBeginUser;
  if (!user->imex) {
    PetscInt           row[] = {0, 1}, col[] = {0};
    PetscScalar        J[2][1];
    const PetscScalar *u;
    PetscCall(VecGetArrayRead(U, &u));
    J[0][0] = 0;
    J[1][0] = -((1. - u[0] * u[0]) * u[1] - u[0])/(user->mu * user->mu);;
    PetscCall(MatSetValues(A, 2, row, 1, col, &J[0][0], INSERT_VALUES));
    PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
    PetscCall(VecRestoreArrayRead(U, &u));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode MaxStepMonitor(TS ts, PetscInt step, PetscReal time, Vec U, void *ctx)
{
    User user = (User)ctx;
    PetscReal dt;

    PetscFunctionBeginUser;
    PetscCall(TSGetTimeStep(ts, &dt));
    if (dt > user->max_dt) user->max_dt = dt;
    PetscFunctionReturn(0);
}


int main(int argc, char **argv)
{
  TS             ts;
  PetscBool      implicitform = PETSC_FALSE;
  PetscScalar   *x_ptr, *y_ptr, derp;
  PetscMPIInt    size;
  struct _n_User user;
  PetscReal tol = 1.e-7; // default value

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize program
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, help));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCheck(size == 1, PETSC_COMM_WORLD, PETSC_ERR_WRONG_MPI_SIZE, "This is a uniprocessor example only!");

  PetscCall(PetscLogDefaultBegin());
  PetscLogEvent USER_EVENT;
  PetscLogEventRegister("TimeForwardAndAdjoint", TS_CLASSID, &USER_EVENT);
  

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Set runtime options
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  user.next_output = 0.0;
  user.mu          = 0.001;
  user.steps       = 0;
  user.ftime       = 0.5;
  user.imex        = PETSC_FALSE;
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-mu", &user.mu, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-implicitform", &implicitform, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-imexform", &user.imex, NULL));

  // Parse user command-line argument
  PetscOptionsGetReal(NULL, NULL, "-tol", &tol, NULL);

  // Show run type ( either implicit or explicit or imex)
  if (user.imex) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Running in IMEX form\n"));
  } else if (implicitform) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Running in implicit form\n"));
  } else {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Running in explicit form\n"));
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Create necessary matrix and vectors, solve same ODE on every process
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(MatCreate(PETSC_COMM_WORLD, &user.A));
  PetscCall(MatSetSizes(user.A, PETSC_DECIDE, PETSC_DECIDE, 2, 2));
  PetscCall(MatSetFromOptions(user.A));
  PetscCall(MatSetUp(user.A));
  PetscCall(MatCreateVecs(user.A, &user.U, NULL));
  PetscCall(MatDuplicate(user.A, MAT_DO_NOT_COPY_VALUES, &user.B));
  PetscCall(MatCreate(PETSC_COMM_WORLD, &user.Jacp));
  PetscCall(MatSetSizes(user.Jacp, PETSC_DECIDE, PETSC_DECIDE, 2, 1));
  PetscCall(MatSetFromOptions(user.Jacp));
  PetscCall(MatSetUp(user.Jacp));
  PetscCall(MatDuplicate(user.Jacp, MAT_DO_NOT_COPY_VALUES, &user.Jacprhs));
  PetscCall(MatZeroEntries(user.Jacprhs));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create timestepping solver context
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(TSCreate(PETSC_COMM_WORLD, &ts));
  PetscCall(TSSetEquationType(ts, TS_EQ_ODE_EXPLICIT)); /* less Jacobian evaluations when adjoint BEuler is used, otherwise no effect */
  
  if (user.imex) {
    PetscCall(TSSetTolerances(ts, tol, NULL, tol, NULL));
    PetscCall(TSSetIFunction(ts, NULL, IFunction, &user));
    PetscCall(TSSetIJacobian(ts, user.A, user.A, IJacobian, &user));
    PetscCall(TSSetIJacobianP(ts, user.Jacp, IJacobianP, &user));
    PetscCall(TSSetRHSFunction(ts, NULL, RHSFunction, &user));
    PetscCall(TSSetRHSJacobian(ts, user.B, NULL, RHSJacobian, &user));
    PetscCall(TSSetRHSJacobianP(ts, user.Jacprhs, NULL, &user));
    PetscCall(TSSetType(ts, TSARKIMEX));
  } else {
    if (implicitform) {
      PetscCall(TSSetIFunction(ts, NULL, IFunction, &user));
      PetscCall(TSSetIJacobian(ts, user.A, user.A, IJacobian, &user));
      PetscCall(TSSetIJacobianP(ts, user.Jacp, IJacobianP, &user));
      PetscCall(TSSetType(ts, TSCN));
    } else {
      TSSetTolerances(ts, tol, NULL, tol, NULL);
      PetscCall(TSSetRHSFunction(ts, NULL, RHSFunction, &user));
      PetscCall(TSSetRHSJacobian(ts, user.A, user.A, RHSJacobian, &user));
      PetscCall(TSSetRHSJacobianP(ts, user.Jacp, RHSJacobianP, &user));
      PetscCall(TSSetType(ts, TSRK));
      PetscCall(TSRKSetType(ts, TSRK5DP));
    }
  }
  PetscCall(TSSetMaxTime(ts, user.ftime));
  PetscCall(TSSetTimeStep(ts, 0.001));
  PetscCall(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set initial conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecGetArray(user.U, &x_ptr));
  x_ptr[0] = 2.0;
  x_ptr[1] = -2.0 / 3.0 + 10.0 / (81.0 /user.mu) - 292.0 / (2187.0 /(user.mu * user.mu));
  PetscCall(VecRestoreArray(user.U, &x_ptr));
  PetscCall(TSSetTimeStep(ts, 0.001));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Save trajectory of solution so that TSAdjointSolve() may be used
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(TSSetSaveTrajectory(ts));

  // Set trajectory type to basic
  TSTrajectory tj;
  TSGetTrajectory(ts, &tj);
  TSTrajectorySetSolutionOnly(tj, PETSC_TRUE);
  PetscCall(TSTrajectorySetType(tj, ts, TSTRAJECTORYMEMORY));


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Initialuze the adjoint variables
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(MatCreateVecs(user.A, &user.lambda[0], NULL));
  /* Set initial conditions for the adjoint integration */
  PetscCall(VecGetArray(user.lambda[0], &y_ptr));
  y_ptr[0] = 1.0;
  y_ptr[1] = 0.0;
  PetscCall(VecRestoreArray(user.lambda[0], &y_ptr));
  PetscCall(MatCreateVecs(user.A, &user.lambda[1], NULL));
  PetscCall(VecGetArray(user.lambda[1], &y_ptr));
  y_ptr[0] = 0.0;
  y_ptr[1] = 1.0;
  PetscCall(VecRestoreArray(user.lambda[1], &y_ptr));

  PetscCall(MatCreateVecs(user.Jacp, &user.mup[0], NULL));
  PetscCall(VecGetArray(user.mup[0], &x_ptr));
  x_ptr[0] = 0.0;
  PetscCall(VecRestoreArray(user.mup[0], &x_ptr));
  PetscCall(MatCreateVecs(user.Jacp, &user.mup[1], NULL));
  PetscCall(VecGetArray(user.mup[1], &x_ptr));
  x_ptr[0] = 0.0;
  PetscCall(VecRestoreArray(user.mup[1], &x_ptr));
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Store max step size and solution
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  // To store the solution 
  Vec UCopy;
  PetscCall(VecDuplicate(user.U, &UCopy));

  PetscCall(TSSetFromOptions(ts));

  user.max_dt = 0.0;

  PetscCall(TSMonitorSet(ts, MaxStepMonitor, &user, NULL));
  PetscCall(TSSolve(ts, user.U));

  //! Start timing here
  PetscLogEventActivate(USER_EVENT);
  PetscLogEventBegin(USER_EVENT, 0, 0, 0, 0);

  PetscCall(TSSolve(ts, user.U)); 

  PetscCall(VecCopy(user.U, UCopy));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Adjoint model starts here
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscCall(TSSetCostGradients(ts, 2, user.lambda, user.mup));

  PetscCall(TSAdjointSolve(ts));

  //! Timing ends here
  PetscLogEventEnd(USER_EVENT, 0, 0, 0, 0);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Display results
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  const PetscScalar *u_ptr,*mup0_ptr, *mup1_ptr;
  
  /* get the forward solution U at final time */
  PetscCall(VecGetArrayRead(UCopy, &u_ptr));


  /* get the adjoint sensitivities wrt mu */
  PetscCall(VecGetArrayRead(user.mup[0], &mup0_ptr));
  PetscCall(VecGetArrayRead(user.mup[1], &mup1_ptr));

  /* print: tolerance U[0] U[1] mup0[0] mup1[0] all space‐separated */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Data is: %.7e %.15e %.15e %.15e %.15e\n",
  tol,
  (double)PetscRealPart(u_ptr[0]),
  (double)PetscRealPart(u_ptr[1]),
  (double)PetscRealPart(mup0_ptr[0]),
  (double)PetscRealPart(mup1_ptr[0])));

  PetscEventPerfInfo eventInfo;
  const char        *name;
  PetscCall(PetscLogEventGetName(USER_EVENT, &name));
  // Display the name of the event
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Event name: %s\n", name));

  // Get the event information
  PetscLogEventInfo event_info;
  PetscCall(PetscLogEventGetPerfInfo(0, USER_EVENT, &eventInfo));

  // Print eventInfo.time with several decimal places
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Time taken for the event: %.15f seconds\n", (double)eventInfo.time));
  PetscPrintf(PETSC_COMM_WORLD, "Max time step used: %g\n", (double)user.max_dt);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(MatDestroy(&user.A));
  PetscCall(MatDestroy(&user.B));
  PetscCall(MatDestroy(&user.Jacp));
  PetscCall(MatDestroy(&user.Jacprhs));
  PetscCall(VecDestroy(&user.U));
  PetscCall(VecDestroy(&user.lambda[0]));
  PetscCall(VecDestroy(&user.lambda[1]));
  PetscCall(VecDestroy(&user.mup[0]));
  PetscCall(VecDestroy(&user.mup[1]));
  PetscCall(TSDestroy(&ts));

  PetscCall(PetscFinalize());
  return 0;
}
