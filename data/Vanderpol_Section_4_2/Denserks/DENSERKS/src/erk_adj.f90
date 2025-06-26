!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Copyright (c) 2007, Mihai Alexe, Adrian Sandu
!                     Virginia Polytechnic Institute and State University
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE RKADJ_MODULE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Module containing the first order adjoint model integrator(subroutine RKINT_ADJ)
!~~~> The user should not call RKINT_ADJ directly. Instead, he/she should call the wrapper
! subroutine RKINT_ADJDR that automatically handles any necessary checkpoint restores
! and forward trajectory recomputations.
!~~~> The checkpointing mechanism is completely transparent to RKINT_ADJ. No checkpoints
! are read, all the forward trajectory information is assumed to be in the memory buffers.
! This allows RKINT_ADJDR (in ./src/erk_adjdr.f90) to handle all checkpoint reads and
! recomputations.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~> use the memory buffers (tapes) implementation
  USE TAPES
!~~~> use the BLAS-like utility functions for vector computations
  USE UTILS

!~~~> use the coefficients for the Runge-Kutta methods defined in RKINT (avoids code duplication)
!~~~> maximum number of Runge-Kutta stages = 16 (for DOPRI8(6) with dense output)
!~~~> A (16 by 16 array, double precision),
!     b (16 by 1 array, double precision),
!     c (16 by 1 array, double precision) : coefficients for the Runge-Kutta methods
!~~~> e (16 by 1 array, double precision) = error coefficients
!~~~> Runge-Kutta Butcher tableau :
!           c | A
!           __|___
!             | b
  USE RK_MODULE, ONLY : s_Max,A,b,c,d,e
  IMPLICIT NONE
  PUBLIC
  SAVE
!~~~>  Statistics on the work performed by the RK method
  INTEGER, PARAMETER :: adNrhs=1, & !ISTATUS(adNrhs) = number of first order adjoint ODE right-hand side evaluations
                        adNstp=2, & !ISTATUS(adNstp) = total number of time steps taken
			adNacc=3, & !ISTATUS(adNacc) = total number of accepted time steps
                        adNrej=4, & !ISTATUS(adNrej) = total number of rejected time steps
			adNhess=5, & !ISTATUS(adNhess) = total number of Hessian-vector products (for 5th order Hermite)
			adNrhsp = 6, & !ISTATUS(adNrhsp) = number of quadrature ODE right-hand side evaluations
                        adNtexit=1, & !RSTATUS(Ntexit) = Texit, the time corresponding to the
                                    !computed forward model solution upon return
			adNhexit=2, & !RSTATUS(adNhexit) = Hexit, last accepted step before exit
			adNhnew=3 !RSTATUS(adNhnew) = Hnew, last predicted step (not yet taken)

 CONTAINS


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE RKINT_ADJ(N, ADY, ADJ_RHS, JACV,          &
           P, W, QUAD_RHS,                           &
           Tstart, Tend,                             &
           AbsTol, RelTol,                           &
           RCNTRL, ICNTRL, RSTATUS, ISTATUS, IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Given an ODE y' = F(y,t), Tstart <= t <= Tend, RKADJ integrates the
!    corresponding continuous adjoint final-value problem:
!                      lambda' = -(dF/dy)*lambda - dg/dy
!                                         Tend >= t >= Tstart
!                      w = -(dF/dp)*lambda - dg/dp
!    using an explicit Runge-Kutta method. Either dense output or Hermite
!    interpolation is used when interpolating the forward solution y(t)
!    during the adjoint model run.
!
!    For details on Runge-Kutta methods and their implementation consult:
!      E. Hairer and G. Wanner
!      "Solving ODEs I. Non-Stiff Problems".
!      Springer series in computational mathematics, Springer-Verlag, 1996.
!    The codes contained in the book inspired this implementation.
!
!    Contact: Mihai Alexe mihai@vt.edu
!             Adrian Sandu asandu@cs.vt.edu
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>   INPUT ARGUMENTS:
!
!      N (integer)      = size of the adjoint ODE system
!
!      ADY(N by 1 array, double precision) = vector of adjoint final conditions (at T=Tend)
!
!      ADJ_RHS (external) = name (external) of the subroutine computing the RHS of
!                the adjoint ODE for a given adjoint state value. This subroutine must
!                have the following form:
!                     SUBROUTINE ADJ_RHS(N, T1, Y1, ADY1, RHS1)
!                       **Inputs:**
!                       INTEGER :: N !costate (adjoint state) size
!                       DOUBLE PRECISION :: T1 !current time moment
!                       DOUBLE PRECISION :: Y1(N) !forward model state at time T1
!                       DOUBLE PRECISION :: ADY1(N) !costate at time T1
!                       **Output:**
!                       DOUBLE PRECISION :: RHS1(N) !RHS of the adjoint ODE
!
!      JACV (external) = name (external) of subroutine computing the value of the Jacobian of
!            F times a given user vector v, i.e. Jac*v. This is needed for the 5th
!            order Hermite interpolation in the adjoint computation. In this case,
!            the procedure is used to compute an approximation of the second
!            derivative of the forward solution: y''(t). If 3rd order
!            Hermite interpolation or dense output are used, the user can pass a
!            dummy pointer for this argument. Otherwise, the subroutine myst have the
!            following form:
!                     SUBROUTINE JACV(N, T1, Y1, F1, JACV1)
!                       **Inputs:**
!                       INTEGER :: N !forward model state size
!                       DOUBLE PRECISION :: T1 !current time moment
!                       DOUBLE PRECISION :: Y1(N), F1(N) !forward model state and
!                                                        !forward ODE RHS at time T1
!                       **Output:**
!                       DOUBLE PRECISION :: JACV1(N)  !Jacobian (evaluated at Y1) times F1
!
!      QUAD_RHS (external) = name(external) of the subroutine computing the RHS of the quadrature
!                 ODE (if present). If the user does not wish to compute sensitivities w.r.t
!                 any parameters, then he/she can pass in a dummy pointer for this argument. 
!                 Otherwise, the subroutine must have the following form:
!                     SUBROUTINE QUAD_RHS(N, P, T1, Y1, ADY1, W1)
!                       **Inputs:**
!                       INTEGER :: N,P !costate size and number of parameters (if any)
!                       DOUBLE PRECISION :: T1 !current time moment
!                       DOUBLE PRECISION :: Y1(N) !forward model state at T1
!                       DOUBLE PRECISION :: ADY1(N) !adjoint model state at T1
!                       **Output:**
!                       DOUBLE PRECISION :: W1(P) !RHS of the quadrature ODE
!
!             The user is advised to obtain the code for the evaluation of Jacobian-vector
!             and Jacobian-transpose-vector products by using the forward mode and the reverse
!             mode of automatic differentiation, respectively. This is computationally efficient
!             in that the full Jacobian is not explicitly formed or stored. Also, the resulting
!             procedure does not incur truncation errors (such errors are inherent with any
!             finite difference approximation).
!
!      P (integer)      =  size of the quadrature ODE system (i.e. number of parameters)
!
!      W (P by 1 array, double precision)  =  quadrature ODE system state at (T->Tend)
!
!      Tstart (double precision), Tend (double precision) :
!            [Tstart,Tend]  = time range of integration. RK_ADJ requires that Tstart <= Tend
!
!    RelTol (N+P by 1 array, double precision),
!    AbsTol (N+P by 1 array, double precision) = user precribed tolerances for the first order adjoint model
!                               variables(RelTol(1:N), AbsTol(1:N)) and quadrature variables
!                               (RelTol(N+1:N+P), AbsTol(N+1:N+P))
!
!     ICNTRL(1:20)    = integer input parameters
!     RCNTRL(1:20)    = real input parameters
!     For more information about ICNTRL and RCNTRL, see the "INPUT PARAMETERS" section below.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     INPUT PARAMETERS:
!
!    Note: For input parameters equal to zero the default values of the
!       corresponding variables are used.
!
!    ICNTRL(1) = 0: AbsTol, RelTol are N-dimensional vectors
!              = 1: AbsTol, RelTol are scalars
!
!    ICNTRL(2)  -> Maximum number of adjoint integration steps
!        (For ICNTRL(2)=0) the default value of 100000 is used
!
!    ICNTRL(3) -> Runge-Kutta method that has been used in the forward model
!                integration. This is needed for the Hermite interpolation or
!                dense output routines (to avoid stage recomputations or
!                additional function calls)
!              = 0 :    DOPRI5(4) (default)
!              = 1 :    Fehlberg's RK2(3)
!              = 2 :    RK3(2) (based on Heun's RK3 scheme)
!              = 3 :    RK4(3) (based on the "3/8-rule")
!              = 4 :    DOPRI5(4)
!              = 5 :    RK6(5)9FM
!              = 6 :    DOPRI8(6)
!
!    ICNTRL(4) -> Interpolation method used in the adjoint solve.
!              = 1: 3rd order Hermite
!              = 2: 5th order Hermite
!              else : Dense output (works **only** for Dopri5(4), RK6(5) and Dopri8(6))
!
!    ICNTRL(5) -> The forward ODE is dependent on a set of parameters or parameter-free
!             = 0: no parameters are present
!             = 1: the forward ODE depends on a set of P parameters
!
!    RCNTRL(1)  -> Hmin, lower bound for the integration step size
!          It is strongly recommended to keep Hmin = ZERO
!    RCNTRL(2)  -> Hmax, upper bound for the integration step size
!    RCNTRL(3)  -> Hstart, starting value for the integration step size
!
!    RCNTRL(4)  -> FacMin, lower bound on step decrease factor (default=0.2)
!    RCNTRL(5)  -> FacMax, upper bound on step increase factor (default=5)
!    RCNTRL(6)  -> FacRej, step decrease factor after multiple rejections
!                          (default=0.1)
!    RCNTRL(7)  -> FacSafe, by which the new step is slightly smaller
!         than the predicted value  (default=0.9)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!    OUTPUT ARGUMENTS:
!    -----------------
!
!    T (double precision)    -> T value for which the solution has been computed
!                     (after successful return T=Tstart).
!
!    ADY (N by 1 array, double precision) -> Adjoint model numerical solution at T
!
!    W (P by 1 array, double precision)  -> Quadrature ODE system solution at T
!
!    IERR (integer)    -> Reports on successfulness upon return:
!                    = 1 for success
!                    < 0 for error (value equals error code)
!
!    ISTATUS (20 by 1 array, integer)
!    RSTATUS (20 by 1 array, double precision) :
!
!    ISTATUS(1)  -> No. of adjoint ODE RHS evaluations
!    ISTATUS(2)  -> No. of steps
!    ISTATUS(3)  -> No. of accepted steps
!    ISTATUS(4)  -> No. of rejected steps (except at very beginning)
!    ISTATUS(5)  -> No. of Jacobian times vector products
!    ISTATUS(6)  -> No. of quadrature ODE RHS evaluations (=0 if n/a)
!
!    RSTATUS(1)  -> Texit, the time corresponding to the
!                     computed Y upon return
!    RSTATUS(2)  -> Hexit, last accepted step before exit
!    RSTATUS(3)  -> Hnew, last predicted step (not yet taken)
!                   For multiple restarts, use Hnew as Hstart
!                     in the subsequent run
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

!~~~>  Arguments
   INTEGER,          INTENT(IN) :: N
   INTEGER,          INTENT(IN) :: P
   DOUBLE PRECISION, INTENT(INOUT) :: ADY(N)
   DOUBLE PRECISION, INTENT(INOUT) :: W(P)
   DOUBLE PRECISION, INTENT(INOUT) :: Tstart,Tend
   DOUBLE PRECISION, INTENT(IN)    :: AbsTol(N+P),RelTol(N+P)
   INTEGER,          INTENT(IN)    :: ICNTRL(20)
   DOUBLE PRECISION, INTENT(IN)    :: RCNTRL(20)
   INTEGER,          INTENT(INOUT) :: ISTATUS(20)
   DOUBLE PRECISION, INTENT(INOUT) :: RSTATUS(20)
   INTEGER,          INTENT(OUT)   :: IERR
   EXTERNAL ADJ_RHS, JACV, QUAD_RHS

!~~~>  Local variables
   DOUBLE PRECISION :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   DOUBLE PRECISION :: Hmin, Hmax, Hstart
   DOUBLE PRECISION :: Texit
   INTEGER       :: i, UplimTol, Max_no_steps
   INTEGER       :: FwdMethod, IntMethod
   LOGICAL       :: VectorTol
   LOGICAL       :: Parameters

!~~~>   Parameters
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.d0, ONE  = 1.d0
   !this is used as the initial time step if no time step is provided in RCNTRL(3)
   DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0d-5
   INTEGER, PARAMETER :: RK2=1, RK3=2, RK4=3, RK5=4, RK6=5, RK8=6

   !number of stages
   INTEGER :: rk_S
   !order of the RK method
   DOUBLE PRECISION :: rk_ELO

!~~~>  Initialize statistics
   ISTATUS(1:8) = 0
   RSTATUS(1:3) = ZERO

!~~~>  For Scalar tolerances (ICNTRL(1) .NE. 0)  the code uses AbsTol(1) and RelTol(1)
!~~~>  For Vector tolerances (ICNTRL(1) == 0) the code uses AbsTol(1:N+P) and RelTol(1:N+P)
   IF (ICNTRL(1) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N+P
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(2) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(2) > 0) THEN
      Max_no_steps=ICNTRL(2)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(2)=',ICNTRL(2)
      CALL rk_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~> Determine the RK method which will be used
   IF ((ICNTRL(3) .EQ. 0) .OR. (ICNTRL(3) .EQ. 4)) THEN
      !Dopri5(4) with 4th order dense output
      FwdMethod = RK5
      rk_S = 7
      rk_ELO = 5.d0
   ELSEIF (ICNTRL(3) .EQ. 1) THEN
      !RK2(3)
      FwdMethod = RK2
      rk_S = 3
      rk_ELO = 2.d0
   ELSEIF (ICNTRL(3) .EQ. 2) THEN
      !RK3(2)
      FwdMethod = RK3
      rk_S = 4
      rk_ELO = 3.d0
   ELSEIF (ICNTRL(3) .EQ. 3) THEN
      !RK4(3)
      FwdMethod = RK4
      rk_S = 6
      rk_ELO = 4.D0
   ELSEIF (ICNTRL(3) .EQ. 5) THEN
      !RK6(5) w/ 6th order dense output
      FwdMethod = RK6
      rk_S = 12
      rk_ELO = 6.d0
   ELSEIF (ICNTRL(3) .EQ. 6) THEN
      !RK8(6) with 7th order dense output
      FwdMethod = RK8
      rk_S = 16
      rk_ELO = 8.d0
   ELSE
      !Incorrect user selection
      PRINT * ,'Runge-Kutta method not implemented: ICNTRL(3)=',ICNTRL(3)
      CALL rk_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~> Interpolation method
   IF ((FwdMethod .EQ. RK2) .OR. (FwdMethod .EQ. RK3) .OR. &
        &(FwdMethod .EQ. RK4)) THEN
       !For lower-order methods (up to 4), only Hermite interpolation is available
       IF (.NOT. ((ICNTRL(4) .EQ. 1) .OR. (ICNTRL(4) .EQ. 2))) THEN
            PRINT * ,'Wrong choice of interpolation method for the selected ', & 
                   & ' RK method ICNTRL(4)=',ICNTRL(4)
            CALL rk_ErrorMsg(-10,Tstart,ZERO,IERR)
            RETURN
       END IF
   ELSE
       !Higher order methods also support dense output
       IF (.NOT. ((ICNTRL(4) .EQ. 1) .OR. (ICNTRL(4) .EQ. 2) .OR. &
           & (ICNTRL(4) .EQ. 0))) THEN
            PRINT * ,'Wrong choice of interpolation method for the selected ', & 
                   & ' RK method ICNTRL(4)=',ICNTRL(4)
            CALL rk_ErrorMsg(-10,Tstart,ZERO,IERR)
            RETURN
       END IF
   END IF

   IntMethod = ICNTRL(4)

!~~~> Check if the problem is parameter-dependent
   IF (ICNTRL(5) == 0) THEN
      Parameters = .FALSE.
   ELSEIF (ICNTRL(5) == 1) THEN
      Parameters = .TRUE.
   ELSE
      PRINT * ,'Invalid parameters option: ICNTRL(5)=',ICNTRL(5)
      CALL rk_ErrorMsg(-11,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~> Unit roundoff (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~> Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL rk_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL rk_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL rk_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2d0
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL rk_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 5.d0
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL rk_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1d0
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL rk_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9d0
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL rk_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) &
         .OR. (RelTol(i) <= 10.d0*Roundoff) &
         .OR. (RelTol(i) >= 1.0) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL rk_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO

!~~~>   Make sure that Tstart <= Tend; RKINT_ADJ will integrate
! the first order adjoint model from Tend to Tstart
   IF (Tstart .gt. Tend) THEN
       PRINT * ,'User-selected Tstart > Tend.'
       CALL rk_ErrorMsg(-8,Tstart,ZERO,IERR)
       RETURN
   END IF

!~~~>  CALL Runge-Kutta core integrator method
   CALL RKAdj_Integrator(ADY, ADJ_RHS,              &
        W, QUAD_RHS,                                &
        Tstart, Tend, Texit,                        &
        AbsTol, RelTol,                             &
!  Integration parameters
        VectorTol, Max_no_steps,                    &
        Parameters,                                 &
        Roundoff, Hmin, Hmax, Hstart,               &
        FacMin, FacMax, FacRej, FacSafe,            &
!  Real and integer output parameters
        RSTATUS, ISTATUS,                           &
!  Error indicator
        IERR)

 CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE rk_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages. 
! Inputs:
!      Code (integer) : Error code
!      T (double precision) : Time at which the error occured
!      H (double precision) : Last time step taken before the error occured
! Output:
!      IERR (integer) : Error indicator. Upon exit from rk_ErrorMsg, IERR = Code
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR

   IERR = Code
   PRINT * , &
     'Forced exit from RKADJ due to the following error:' 

   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)
      PRINT * , '--> Improper selection of Runge-Kutta method.'
    CASE (-3)
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4) 
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T or H < Roundoff'
    CASE (-8)
      PRINT * , '--> User selected Tstart > Tend. ', &
                'RKADJ requires that Tstart <= Tend.'
    CASE (-9)
      PRINT * , '--> Allocation error! Check available physical memory'
    CASE (-10)
      PRINT * , '--> Improper choice of interpolation strategy.'
    CASE (-11)
      PRINT * , '--> Invalid option for the parameters (must be either 0 or 1).'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT

   PRINT *, "T=", T, "and H=", H

 END SUBROUTINE rk_ErrorMsg



SUBROUTINE RKAdj_Integrator (ADY, ADJ_RHS,              &
        W, QUAD_RHS,                                    &
        Tstart, Tend, T,                                &
        AbsTol, RelTol,                                 &
!~~~> Integration parameters
        VectorTol, Max_no_steps, Parameters,            &
        Roundoff, Hmin, Hmax, Hstart,                   &
        FacMin, FacMax, FacRej, FacSafe,                &
        RSTATUS, ISTATUS,                               &
!~~~> Error indicator
        IERR )

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Core Runge-Kutta integrator. Performs the actual integration of the 
! 1st order adjoint ODE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Input parameters:
!         ADY (N by 1 array, double precision) :
!                    first order adjoint model state at Tend
!         W (P by 1 array, double precision) :
!                    quadrature ODE system state at Tend
!         ADJ_RHS (external) : pointer to (user-supplied) subroutine that evaluates
!                          the right-hand side of the first order adjoint model ODE
!         QUAD_RHS (external) : pointer to (user-supplied) subroutine that evaluates
!                          the right-hand side of the quadrature ODE system
!         Tend (double precision) : integration start time
!         Tstart (double precision) : integration end time
!         AbsTol (N+P by 1 array, double precision) : absolute tolerances
!         RelTol (N+P by 1 array, double precision) : relative tolerances
!         VectorTol (logical) : vector tolerances (TRUE) or scalar tolerances (FALSE)
!         Max_no_steps (integer) : maximum number of integration steps
!         Parameters (logical) : Flag signaling whether or not the system is parameter-dependent
!                       = TRUE : System is dependent on P parameters
!                       = FALSE : No parameters present
!                       (this allows the integrator to avoid parameter-related computations when
!                        parameters are not present)
!         Roundoff (double precision) : machine epsilon
!         Hmin (double precision) : lower bound for the integration step size
!         Hmax (double precision) : upper bound for the integration step size
!         Hstart (double precision) : starting value for the integration step size
!         FacMin (double precision) : lower bound on step decrease factor
!         FacMax (double precision) : upper bound on step increase factor
!         FacRej (double precision) : step decrease factor after multiple rejections
!         FacSafe (double precision) : factor by which the new step is slightly smaller
!                                      than the predicted value
!
!
! Output parameters:
!         T (double precision) : time value for which the solution has been computed
!                     (after successful return T=Tstart) 
!         ADY (N by 1 array, double precision) :
!                    first order adjoint model state at T
!         W (P by 1 array, double precision) :
!                    quadrature ODE system state at T
!         Nc (integer) : total number of file checkpoints written during the integration
!         RSTATUS (20 by 1 array, double precision) : integrator output parameters
!         ISTATUS (20 by 1 array, double precision) : integrator output parameters (statistics)
!
!         For more information on RSTATUS and ISTATUS, see the comments describing the RKINT subroutine
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



   IMPLICIT NONE
   EXTERNAL ADJ_RHS
   EXTERNAL QUAD_RHS
!~~~> Input: the initial condition at Tstart; Output: the solution at T
   DOUBLE PRECISION, INTENT(INOUT) :: ADY(N)
!~~~> Input: the gradient of the cost function w.r.t. the parameters
   DOUBLE PRECISION, INTENT(INOUT) :: W(P)
!~~~> Input: integration interval
   DOUBLE PRECISION, INTENT(INOUT) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   DOUBLE PRECISION, INTENT(OUT) ::  T
!~~~> Input: tolerances
   DOUBLE PRECISION, INTENT(IN) ::  AbsTol(N+P), RelTol(N+P)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: VectorTol
   LOGICAL, INTENT(IN) :: Parameters
   DOUBLE PRECISION, INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   DOUBLE PRECISION, INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   DOUBLE PRECISION, INTENT(INOUT) :: RSTATUS(20)

!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR

! ~~~~ Local variables
   !New adjoint solution
   DOUBLE PRECISION :: ADYnew(N)
   !New gradient
   DOUBLE PRECISION :: Wnew(P)
   !Kadj = the stages for the integration of the adjoint equations
   DOUBLE PRECISION :: Kadj(N,rk_S)
   !Kw = the stages for the quadrature equations
   DOUBLE PRECISION :: Kw(P,rk_S)
   !Error vectors for the quadratures and adjoint variables
   DOUBLE PRECISION :: Werr(P), ADYerr(N)
   !Enlarged error vectors when the quadrature equations are added
   !This allows full error control (quadrature variables participate in the error control mechanism)
   DOUBLE PRECISION :: ADYW(N+P), ADYWnew(N+P), ADYWerr(N+P)

   DOUBLE PRECISION :: H, Hnew, Fac, Tau
   DOUBLE PRECISION :: Err, Yval(N)
   DOUBLE PRECISION :: Jac0(N), Jac(N), JacP0(P), JacP(P)
   INTEGER :: j, istage
   INTEGER :: allocErr = 0
   LOGICAL :: RejectLastH, RejectMoreH
!~~~>  Local parameters (constants)
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.d0, ONE  = 1.d0
   DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0d-5

!~~~> Check if the tapes have been allocated (the user has to take care of this in the driver
! by calling rk_AllocateTapes )
   CALL tapes_CheckAllocated(allocErr)
   IF (allocErr .NE. 1) THEN
        CALL rk_ErrorMsg(-2,Tend,ZERO,IERR)
   END IF

!~~~>  Initial preparations for the adjoint integration
   T = Tend
   RSTATUS(adNhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.d0*Roundoff) H = DeltaMin

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

AdjTimeLoop: DO WHILE ((Tstart-T)+Roundoff <= ZERO)

   IF ( ISTATUS(adNstp) > Max_no_steps ) THEN  ! Too many steps
      CALL rk_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1d0*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL rk_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tstart
   H = MIN(H,ABS(T-Tstart))

!~~~>  Repeat step calculation until current step accepted
AdjUntilAccepted: DO

   CALL rk_InterpolateFwdSol(T, Yval)

   CALL ADJ_RHS(N,T,Yval,ADY,Jac0)
   !scale the RHS by -1 because of negative timesteps
   CALL WSCAL(N,-1.d0,Jac0,1)
   ISTATUS(adNrhs) = ISTATUS(adNrhs) + 1

   IF (Parameters) THEN
      CALL QUAD_RHS(N,P,T,Yval,ADY,JacP0)
      !scale the RHS by -1 because of negative timesteps
      CALL WSCAL(P,-1.d0,JacP0,1)
      ISTATUS(adNrhsp) = ISTATUS(adNrhsp) + 1
   END IF

!~~~>   Compute the stages
AdjStage: DO istage = 1, rk_S

     IF (istage .EQ. 1) THEN

         !Jac = Jac0
         CALL WCOPY(N,Jac0,1,Jac,1)
         IF (Parameters) THEN
            !JacP = JacP0
            CALL WCOPY(P,JacP0,1,JacP,1)
         END IF

     ELSE

         !ADYnew = ADY
         CALL WCOPY(N,ADY,1,ADYnew,1)

         IF (Parameters) THEN
            !Wnew = W
            CALL WCOPY(P,W,1,Wnew,1)
         END IF

         DO j = 1, istage-1
             !ADYnew = ADYnew + H*A(istage,j)*Kadj(:,j)
             CALL WAXPY(N,H*A(istage,j),Kadj(:,j),1,ADYnew,1)

             IF (Parameters) THEN
                 !Wnew = Wnew + H*A(istage,j)*Kw(:,j)
                 CALL WAXPY(P,H*A(istage,j),Kw(:,j),1,Wnew,1)
             END IF

         END DO

         Tau = T - c(istage)*H

         CALL rk_InterpolateFwdSol(Tau,Yval)

         CALL ADJ_RHS(N,Tau,Yval,ADYnew,Jac)
         CALL WSCAL(N,-1.d0,Jac,1)
         ISTATUS(adNrhs) = ISTATUS(adNrhs) + 1

         IF (Parameters) THEN
              CALL QUAD_RHS(N,P,Tau,Yval,ADYnew,JacP)
              CALL WSCAL(P,-1.d0,JacP,1)
              ISTATUS(adNrhsp) = ISTATUS(adNrhsp) + 1
         END IF

     END IF

     !Kadj(:,istage) = Jac
     CALL WCOPY(N,Jac,1,Kadj(:,istage),1)

     IF (Parameters) THEN
         !Kw(:,istage) = JacP
         CALL WCOPY(P,JacP,1,Kw(:,istage),1)
     END IF

   END DO AdjStage


!~~~>  Compute the new solution
   !ADYnew = ADY
   CALL WCOPY(N,ADY,1,ADYnew,1)
   DO j=1,rk_S
        !ADYnew = ADYnew + H*b(j)*Kadj(:,j)
        CALL WAXPY(N,H*b(j),Kadj(:,j),1,ADYnew,1)
   END DO

   IF (Parameters) THEN
!~~~> Compute the new gradient Wnew
        !Wnew = W
        CALL WCOPY(P,W,1,Wnew,1)
        DO j=1,rk_S
            !Wnew = Wnew + H*b(j)*Kw(:,j)
            CALL WAXPY(P,H*b(j),Kw(:,j),1,Wnew,1)
        END DO
   END IF

!~~~>  Compute the error estimation
   !ADYerr = 0.d0
   CALL SET2ZERO(N,ADYerr)
   DO j=1,rk_S
        !ADYerr = ADYerr + H*e(j)*Kadj(:,j)
        CALL WAXPY(N,H*e(j),Kadj(:,j),1,ADYerr,1)
   END DO

   IF (Parameters) THEN 
!~~~>  Compute the error estimation for the quadrature equations solution W
      !Werr = 0.d0
      CALL SET2ZERO(P,Werr)
      DO j=1,rk_S
          !Werr = Werr + H*e(j)*Kw(:,j)
          CALL WAXPY(P,H*e(j),Kw(:,j),1,Werr,1)
      END DO
   END IF

!~~~> Compute the scaled error norm of ADYerr (no parameters) or ADYWerr (w/ P parameters)
   IF (.NOT. Parameters) THEN 
        Err = rk_ErrorNorm ( N, ADY, ADYnew, ADYerr, AbsTol, RelTol, VectorTol )
   ELSE
        CALL SET2ZERO(N+P,ADYWerr)
        CALL WCOPY(N,ADY,1,ADYW(1:N),1)
        CALL WCOPY(P,W,1,ADYW(N+1:N+P),1)
        CALL WCOPY(N,ADYnew,1,ADYWnew(1:N),1)
        CALL WCOPY(P,Wnew,1,ADYWnew(N+1:N+P),1)
        CALL WCOPY(N,ADYerr,1,ADYWerr(1:N),1)
        CALL WCOPY(P,Werr,1,ADYWerr(N+1:N+P),1)
        Err = rk_ErrorNorm ( N+P, ADYW, ADYWnew, ADYWerr, AbsTol, RelTol, VectorTol )
   END IF

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/rk_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(adNstp) = ISTATUS(adNstp) + 1

   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      !increment the "accepted steps" counter
      ISTATUS(adNacc) = ISTATUS(adNacc) + 1
      !update the solution
      !ADY = ADYnew
      CALL WCOPY(N,ADYnew,1,ADY,1)
      !update the gradient w.r.t. parameters if needed
      IF (Parameters) THEN
          !W = Wnew
          CALL WCOPY(P,Wnew,1,W,1)
      END IF

      !update current time
      T = T - H
      !make sure Hmin <= Hnew <= Hmax
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTATUS(adNhexit) = H
      RSTATUS(adNhnew)  = Hnew
      RSTATUS(adNtexit) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      !final update update step size
      H = Hnew
      EXIT AdjUntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED

   ELSE           !~~~> Reject step

      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(adNacc) >= 1)  ISTATUS(adNrej) = ISTATUS(adNrej) + 1

   END IF ! Err <= 1

   END DO AdjUntilAccepted

   END DO AdjTimeLoop

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

END SUBROUTINE RKAdj_Integrator



DOUBLE PRECISION FUNCTION rk_ErrorNorm ( N, Y, Ynew, Yerr, &
                               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
! Inputs:
!      N (integer) : dimension of the vectors Y, Ynew, Yerr, AbsTol, RelTol
!      Ynew (N by 1 array, double precision)
!             Solution vector
!      Yerr (N by 1 array, double precision)
!             Error vector
!      AbsTol (N by 1 array, double precision)
!             Integrator absolute tolerances
!      RelTol (N by 1 array, double precision)
!             Integrator absolute tolerances
!      VectorTol (logical)
!             = FALSE : Scalar tolerances (only AbsTol(1) and RelTol(1) are
!                                         used in the error norm computations)
!             = TRUE : Vector tolerances
! Outputs: 
!      rk_ErrorNorm (double precision) : the scaled error norm of the vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

! Input arguments
   INTEGER, INTENT(IN) :: N
   DOUBLE PRECISION, INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   DOUBLE PRECISION :: Err, Scale, Ymax
   INTEGER  :: i
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0

   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   rk_ErrorNorm = MAX(Err,1.0d-10)

END FUNCTION rk_ErrorNorm



SUBROUTINE rk_DenseOutput4(T, Yval)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Provides 4th order dense output for the DOPRI5(4) Runge-Kutta method.
!
! Implements the formulas given by Dormand and Prince ("Runge-Kutta Triples",
! Comp & Maths with Applic., Vol 12A, p1007-1017, 1986)
! No additional function evaluations are required. The procedure results
! (globally) in a C1 (continuously differentiable) approximation to the forward
! model solution.
!
! Inputs:
!        T (double precision) : time at which a forward solution approximation
!                               is required
!
! Outputs:
!        Yval (N by 1 array, double precision) : interpolated forward solution at T
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: T

  DOUBLE PRECISION :: Yp(N), Kp(N,s_Max)
  DOUBLE PRECISION :: Tp, Hp, Theta

  DOUBLE PRECISION, INTENT(OUT) :: Yval(N)

  DOUBLE PRECISION, SAVE :: b_theta(7)
  DOUBLE PRECISION :: Temp1, Temp2
  INTEGER :: i

  CALL tapes_LookupDense(T,Tp,Hp,Kp,Yp)

  Theta = (T - Tp) / Hp
  Temp1 = Theta**2.d0 * (Theta-1.d0)**2.d0
  Temp2 = Theta**2.d0 * (3.d0 - 2.d0*Theta)

  !we use the dense output formulas from Hairer's text
  b_theta(1) = Temp2*b(1) + Theta*(Theta-1.d0)**2.d0 - 5.d0*Temp1*(2558722523.d0 - 31403016.d0*Theta)/11282082432.d0
  b_theta(2) = 0.d0
  b_theta(3) = Temp2*b(3) + 100.d0*Temp1*(882725551.d0 - 15701508.d0*Theta)/32700410799.d0
  b_theta(4) = Temp2*b(4) - 25.d0*Temp1*(443332067.d0 - 31403016.d0*Theta)/1880347072.d0
  b_theta(5) = Temp2*b(5) + 32805.d0*Temp1*(23143187.d0 - 3489224.d0*Theta)/199316789632.d0
  b_theta(6) = Temp2*b(6) - 55.d0*Temp1*(29972135.d0 - 7076736.d0*Theta)/822651844.d0
  b_theta(7) = Theta**2.d0 * (Theta-1.d0) + 10.d0*Temp1*(7414447.d0 - 829305.d0*Theta)/29380423.d0

  !Yval = Yc
  CALL WCOPY(N,Yp,1,Yval,1)
  DO i=1,rk_S
      !Yval = Yval + Hc*b_Theta(i)*Kc(:,i)
      CALL WAXPY(N,Hp*b_theta(i),Kp(:,i),1,Yval,1)
  END DO

END SUBROUTINE rk_DenseOutput4



SUBROUTINE rk_DenseOutput6(T, Yval)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Provides 6th order dense output for the RK6(5) integrator
!
! Implements the formulas given by Baker, Dormand, Gilmore and Prince
! ("Continuous approximation with embedded Runge-Kutta methods",Appl. Num. Math. 22, 1996)
!
! Inputs:
!        T (double precision) : time at which a forward model solution approximation
!                               is required
!
! Outputs:
!        Yval (N by 1 array, double precision) : interpolated forward model solution at T
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: T

  DOUBLE PRECISION :: Yp(N), Kp(N,s_Max)
  DOUBLE PRECISION :: Tp, Hp, Theta1

  DOUBLE PRECISION, INTENT(OUT) :: Yval(N)

  DOUBLE PRECISION, SAVE :: b_theta(12,6), bb(12)
  DOUBLE PRECISION, SAVE :: theta(6)
  INTEGER :: i


  CALL tapes_LookupDense(T,Tp,Hp,Kp,Yp)
  Theta1 = (T - Tp) / Hp

  b_theta = 0.d0
  b_theta(1,1) = -2071.d0/300.d0
  b_theta(4,1) = -483328.d0/23595.d0
  b_theta(5,1) = -531441.d0/10285.d0
  b_theta(6,1) = 8576.d0/235.d0
  b_theta(7,1) = -1977326743.d0/75409620.d0
  b_theta(8,1) = 259.d0/15.d0
  b_theta(9,1) = 1112.d0/175.d0
  b_theta(10,1) = 448.d0/5.d0
  b_theta(11,1) = 864.d0/25.d0
  b_theta(12,1) = -13824.d0/175.d0

  b_theta(1,2) = 13811.d0/600.d0
  b_theta(4,2) = 241664.d0/4719.d0
  b_theta(5,2) = 531441.d0/4114.d0
  b_theta(6,2) = -4288.d0/47.d0
  b_theta(7,2) = 1977326743.d0/30163848.d0
  b_theta(8,2) = -259.d0/6.d0
  b_theta(9,2) = -2636.d0/175.d0
  b_theta(10,2) = -1264.d0/5.d0
  b_theta(11,2) = -2592.d0/25.d0
  b_theta(12,2) = 41472.d0/175.d0

  b_theta(1,3) = -53923.d0/1800.d0
  b_theta(4,3) = -241664.d0/5445.d0
  b_theta(5,3) = -2302911.d0/20570.d0
  b_theta(6,3) = 55744.d0/705.d0
  b_theta(7,3) = -1977326743.d0/34804440.d0
  b_theta(8,3) = 3367.d0/90.d0
  b_theta(9,3) = 949.d0/75.d0
  b_theta(10,3) = 3767.d0/15.d0
  b_theta(11,3) = 2907.d0/25.d0
  b_theta(12,3) = -6336.d0/25.d0

  b_theta(1,4) = 139189.d0/7200.d0
  b_theta(4,4) = 1147904.d0/70785.d0
  b_theta(5,4) = 3365793.d0/82280.d0
  b_theta(6,4) = -20368.d0/705.d0
  b_theta(7,4) = 37569208117.d0/1809830880.d0
  b_theta(8,4) = -4921.d0/360.d0
  b_theta(9,4) = -2381.d0/525.d0
  b_theta(10,4) = -1534.d0/15.d0
  b_theta(11,4) = -1494.d0/25.d0
  b_theta(12,4) = 19584.d0/175.d0

  b_theta(1,5) = -18487.d0/2880.d0
  b_theta(4,5) = -30208.d0/14157.d0
  b_theta(5,5) = -177147.d0/32912.d0
  b_theta(6,5) = 536.d0/141.d0
  b_theta(7,5) = -1977326743.d0/723932352.d0
  b_theta(8,5) = 259.d0/144.d0
  b_theta(9,5) = 62.d0/105.d0
  b_theta(10,5) = 43.d0/3.d0
  b_theta(11,5) = 63.d0/5.d0
  b_theta(12,5) = -576.d0/35.d0

  b_theta(1,6) = 1.d0

  DO i=1,6
     Theta(i) = Theta1**(6-i)
  END DO

  DO i=1,rk_S
     bb(i) = WDOT(N,b_theta(i,:),1,Theta,1)
  END DO

  CALL WCOPY(N,Yp,1,Yval,1)
  DO i=1,rk_S
      CALL WAXPY(N,Hp*Theta1*bb(i),Kp(:,i),1,Yval,1)
  END DO

END SUBROUTINE rk_DenseOutput6



SUBROUTINE rk_DenseOutput7(T, Yval)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 7th order dense output for the DOPRI8(6) integrator
!
! Implements the dense output following the example of E. Hairer, S.P. Norsett and
! G. Wanner in "Solving ordinary differential equations I: Nonstiff Problems"
! (Springer. 1993).
! Alternatively, the coefficients of the 7th order polynomial b(theta) can be
! obtained by solving a linear system derived from the order conditions. This
! requires that 4 more stages be added to the DOPRI8(6) method.
!
! No additional function or Jacobian evaluations are required. We reuse the
! forward trajectory values stored during the forward model integration.
!
! Inputs:
!        T (double precision) : time at which a forward model solution 
!                               approximation is required
!
! Outputs:
!        Yval (N by 1 array, double precision) : interpolated forward model 
!                                                solution at T
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(INOUT) :: T
  DOUBLE PRECISION, INTENT(OUT) :: Yval(N)

  INTEGER :: j,jj
  DOUBLE PRECISION :: Yc(N), Ycp1(N), Kc(N,rk_S)
  DOUBLE PRECISION :: Theta, Theta1, Tc, Hc
  DOUBLE PRECISION :: Cont(N,0:7)

  !do tape lookup
  CALL tapes_LookupDense(T,Tc,Hc,Kc,Yc,Ycp1)

  !Cont(:,0) = Yc
  CALL WCOPY(N,Yc,1,Cont(:,0),1)

!  Cont(:,1) = Ycp1 - Yc
  CALL WCOPY(N,Yc,1,Cont(:,1),1)
  CALL WSCAL(N,-1.d0,Cont(:,1),1)
  CALL WAXPY(N,1.d0,Ycp1,1,Cont(:,1),1)

!   Cont(:,2) = Hc*Kc(:,1) - Cont(:,1)
  CALL WCOPY(N,Cont(:,1),1,Cont(:,2),1)
  CALL WSCAL(N,-1.d0,Cont(:,2),1)
  CALL WAXPY(N,Hc,Kc(:,1),1,Cont(:,2),1)

!  Cont(:,3) = Cont(:,1) - Hc*Kc(:,13) - Cont(:,2)
  CALL WCOPY(N,Cont(:,2),1,Cont(:,3),1)
  CALL WCOPY(N,Kc(:,13),1,Cont(:,4),1)
  CALL WSCAL(N,-1.d0,Cont(:,3),1)
  CALL WSCAL(N,-1.d0*Hc,Cont(:,4),1)
  CALL WAXPY(N,1.d0,Cont(:,1),1,Cont(:,3),1)
  CALL WAXPY(N,1.d0,Cont(:,4),1,Cont(:,3),1)

  !Cont(:,4:7) = 0.d0
  DO j=4,7
     CALL SET2ZERO(N,Cont(:,j))
  END DO

  DO jj=1,4
     DO j = 1,rk_S
           !Cont(:,3+jj) = Cont(:,3+jj) + Hc*d(3+jj,j)*Kc(:,j)
           CALL WAXPY(N,Hc*d(3+jj,j),Kc(:,j),1,Cont(:,3+jj),1)
     END DO
  END DO

  Theta = (T - Tc) / Hc
  Theta1 = 1.d0 - Theta

  CALL SET2ZERO(N,Yval)
  CALL WCOPY(N,Cont(:,6),1,Yval,1)
  CALL WAXPY(N,Theta,Cont(:,7),1,Yval,1)
  CALL WSCAL(N,Theta1,Yval,1)
  CALL WAXPY(N,1.d0,Cont(:,5),1,Yval,1)
  CALL WSCAL(N,Theta,Yval,1)
  CALL WAXPY(N,1.d0,Cont(:,4),1,Yval,1)
  CALL WSCAL(N,Theta1,Yval,1)
  CALL WAXPY(N,1.d0,Cont(:,3),1,Yval,1)
  CALL WSCAL(N,Theta,Yval,1)
  CALL WAXPY(N,1.d0,Cont(:,2),1,Yval,1)
  CALL WSCAL(N,Theta1,Yval,1)
  CALL WAXPY(N,1.d0,Cont(:,1),1,Yval,1)
  CALL WSCAL(N,Theta,Yval,1)
  CALL WAXPY(N,1.d0,Cont(:,0),1,Yval,1)

END SUBROUTINE rk_DenseOutput7



SUBROUTINE rk_Hermite3(T, Yval)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Template for Hermite interpolation of order 3 on the interval [Ta,Tb]
!  Ta <= T <= Tb
! P = c(1) + c(2)*(x-a) + ... + c(4)*(x-a)^3
! P[a,b] = [Ya,Yb], P'[a,b] = [Ja,Jb]
!
! Input:
!        T (double precision) : time at which a forward solution approximation
!                               is required
!
! Output:
!        Yval (N by 1 array, double precision) : interpolated forward solution at T
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(INOUT) :: T
   DOUBLE PRECISION, INTENT(INOUT) :: Yval(N)

   DOUBLE PRECISION :: Ta, Tb
   DOUBLE PRECISION :: Ya(N), Yb(N), Ja(N), Jb(N)

   DOUBLE PRECISION :: T1, amb(3), C(N,4)
   INTEGER :: i,j
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.d0

   CALL tapes_LookupHermite(T,Ta,Tb,Ya,Yb,Ja,Jb)

   amb(1) = 1.0d0/(Ta - Tb)
   DO i=2,3
     amb(i) = amb(i-1)*amb(1)
   END DO

! c(1) = ya;
   CALL WCOPY(N,Ya,1,C(1,1),1)
! c(2) = ja;
   CALL WCOPY(N,Ja,1,C(1,2),1)
! c(3) = 2/(Ta-Tb)*ja + 1/(Ta-Tb)*jb - 3/(Ta - Tb)^2*ya + 3/(Ta - Tb)^2*yb  ;
   CALL WCOPY(N,Ya,1,C(1,3),1)
   CALL WSCAL(N,-3.0*amb(2),C(1,3),1)
   CALL WAXPY(N,3.0*amb(2),Yb,1,C(1,3),1)
   CALL WAXPY(N,2.0*amb(1),Ja,1,C(1,3),1)
   CALL WAXPY(N,amb(1),Jb,1,C(1,3),1)
! c(4) =  1/(Ta-Tb)^2*ja + 1/(Ta-Tb)^2*jb - 2/(Ta-Tb)^3*ya + 2/(Ta-Tb)^3*yb ;
   CALL WCOPY(N,Ya,1,C(1,4),1)
   CALL WSCAL(N,-2.0*amb(3),C(1,4),1)
   CALL WAXPY(N,2.0*amb(3),Yb,1,C(1,4),1)
   CALL WAXPY(N,amb(2),Ja,1,C(1,4),1)
   CALL WAXPY(N,amb(2),Jb,1,C(1,4),1)

   T1 = T - Ta
   CALL WCOPY(N,C(1,4),1,Yval,1)
   CALL WSCAL(N,T1**3,Yval,1)
   DO j = 3,1,-1
     CALL WAXPY(N,T1**(j-1),C(1,j),1,Yval,1)
   END DO

END SUBROUTINE rk_Hermite3



SUBROUTINE rk_Hermite5(T, Yval)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Template for Hermite interpolation of order 5 on the interval [Ta,Tb]
! Ta <= T <= Tb
! P = c(1) + c(2)*(x-Ta) + ... + c(6)*(x-Ta)^5
! P[Ta,Tb] = [Ya,Yb], P'[Ta,Tb] = [Ja,Jb], P"[Ta,Tb] = [Hessa,Hessb]
!
! Hermite5 requires two additional Jacobian-vector products per function call.
! These are computed via calls to JACV.
!
! Inputs:
!        T (double precision) : time at which a forward solution approximation
!                               is required
!
! Outputs:
!        Yval (N by 1 array, double precision) : interpolated forward solution at T
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   IMPLICIT NONE

   DOUBLE PRECISION, INTENT(IN) :: T
   DOUBLE PRECISION, INTENT(OUT) :: Yval(N)

   DOUBLE PRECISION :: Ya(N), Yb(N)
   DOUBLE PRECISION :: Ja(N), Jb(N), Hessa(N), Hessb(N)
   DOUBLE PRECISION :: Ta, Tb
   DOUBLE PRECISION :: T1, amb(5), C(N,6)
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, HALF = 0.5d0
   INTEGER :: i, j

   CALL tapes_LookupHermite(T,Ta,Tb,Ya,Yb,Ja,Jb)

   amb(1) = 1.0d0 / (Ta-Tb)
   DO i=2,5
     amb(i) = amb(i-1)*amb(1)
   END DO

   CALL JACV(N,Ta,Ya,Ja,Hessa)
   ISTATUS(adNhess) = ISTATUS(adNhess) + 1
   CALL JACV(N,Tb,Yb,Jb,Hessb)
   ISTATUS(adNhess) = ISTATUS(adNhess) + 1

! c(1) = ya;
   CALL WCOPY(N,Ya,1,C(1,1),1)
! c(2) = ja;
   CALL WCOPY(N,Ja,1,C(1,2),1)
! c(3) = ha/2;
   CALL WCOPY(N,Hessa,1,C(1,3),1)
   CALL WSCAL(N,HALF,C(1,3),1)

! c(4) = 10*amb(3)*ya - 10*amb(3)*yb - 6*amb(2)*ja - 4*amb(2)*jb  + 1.5*amb(1)*ha - 0.5*amb(1)*hb ;
   CALL WCOPY(N,Ya,1,C(1,4),1)
   CALL WSCAL(N,10.0*amb(3),C(1,4),1)
   CALL WAXPY(N,-10.0*amb(3),Yb,1,C(1,4),1)
   CALL WAXPY(N,-6.0*amb(2),Ja,1,C(1,4),1)
   CALL WAXPY(N,-4.0*amb(2),Jb,1,C(1,4),1)
   CALL WAXPY(N, 1.5*amb(1),Hessa,1,C(1,4),1)
   CALL WAXPY(N,-0.5*amb(1),Hessb,1,C(1,4),1)

! c(5) =   15*amb(4)*ya - 15*amb(4)*yb - 8.*amb(3)*ja - 7*amb(3)*jb + 1.5*amb(2)*ha - 1*amb(2)*hb ;
   CALL WCOPY(N,Ya,1,C(1,5),1)
   CALL WSCAL(N, 15.0*amb(4),C(1,5),1)
   CALL WAXPY(N,-15.0*amb(4),Yb,1,C(1,5),1)
   CALL WAXPY(N,-8.0*amb(3),Ja,1,C(1,5),1)
   CALL WAXPY(N,-7.0*amb(3),Jb,1,C(1,5),1)
   CALL WAXPY(N,1.5*amb(2),Hessa,1,C(1,5),1)
   CALL WAXPY(N,-amb(2),Hessb,1,C(1,5),1)

! c(6) =   6*amb(5)*ya - 6*amb(5)*yb - 3.*amb(4)*ja - 3.*amb(4)*jb + 0.5*amb(3)*ha -0.5*amb(3)*hb ;
   CALL WCOPY(N,Ya,1,C(1,6),1)
   CALL WSCAL(N, 6.0*amb(5),C(1,6),1)
   CALL WAXPY(N,-6.0*amb(5),Yb,1,C(1,6),1)
   CALL WAXPY(N,-3.0*amb(4),Ja,1,C(1,6),1)
   CALL WAXPY(N,-3.0*amb(4),Jb,1,C(1,6),1)
   CALL WAXPY(N, 0.5*amb(3),Hessa,1,C(1,6),1)
   CALL WAXPY(N,-0.5*amb(3),Hessb,1,C(1,6),1)

   T1 = T - Ta
   CALL WCOPY(N,C(1,6),1,Yval,1)
   DO j = 5,1,-1
     CALL WSCAL(N,T1,Yval,1)
     CALL WAXPY(N,ONE,C(1,j),1,Yval,1)
   END DO

END SUBROUTINE rk_Hermite5


SUBROUTINE rk_InterpolateFwdSol(T,Yval)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This subroutine interpolates the forward solution at a given time T
! rk_InterpolateFwdSol calls the appropriate interpolation routine (either Hermite
! or dense output) according to the user selection of the interpolation method. It
! acts as a front-end for the interpolation mechanism.
!
!  If IntMethod = 1 then call rk_Hermite3
!  If IntMethod = 2 then call rk_Hermite5
!  If IntMethod = 3 then call rk_DenseOuput4/6/7 (depending on FwdMethod)
!
! Input:
!        T (double precision) : time at which a forward solution approximation
!                               is required
!
! Output:
!        Yval (N by 1 array, double precision) : interpolated forward solution at T
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(INOUT) :: T
  DOUBLE PRECISION, INTENT(OUT) :: Yval(N)

  IF ((FwdMethod .EQ. RK2) .OR. &
          & (FwdMethod .EQ. RK3) .OR. &
          & (FwdMethod .EQ. RK4)) THEN

        !Hermite interpolation
        IF (IntMethod .EQ. 1) THEN
             CALL rk_Hermite3(T,Yval)
        ELSE IF (IntMethod .EQ. 2) THEN
             CALL rk_Hermite5(T,Yval)
        END IF

  ELSE IF (FwdMethod .EQ. RK5) THEN

        IF (IntMethod .EQ. 1) THEN
             CALL rk_Hermite3(T,Yval)
        ELSE IF (IntMethod .EQ. 2) THEN
             CALL rk_Hermite5(T,Yval)
        ELSE
             CALL rk_DenseOutput4(T,Yval)
        END IF

  ELSE IF (FwdMethod .EQ. RK6) THEN

        IF (IntMethod .EQ. 1) THEN
             CALL rk_Hermite3(T,Yval)
        ELSE IF (IntMethod .EQ. 2) THEN
             CALL rk_Hermite5(T,Yval)
        ELSE
             CALL rk_DenseOutput6(T,Yval)
        END IF

  ELSE IF (FwdMethod .EQ. RK8) THEN

        IF (IntMethod .EQ. 1) THEN
             CALL rk_Hermite3(T,Yval)
        ELSE IF (IntMethod .EQ. 2) THEN
             CALL rk_Hermite5(T,Yval)
        ELSE
             CALL rk_DenseOutput7(T,Yval)
        END IF

  END IF

END SUBROUTINE rk_InterpolateFwdSol

!END of RKINT_ADJ integrator subroutine
END SUBROUTINE RKINT_ADJ


!END OF RKADJ_MODULE
END MODULE RKADJ_MODULE