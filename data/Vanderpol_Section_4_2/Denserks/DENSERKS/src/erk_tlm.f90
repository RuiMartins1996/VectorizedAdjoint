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

!~~~> Module containing the tangent linear model integrator (subroutine RKINT_TLM)
MODULE RKTLM_MODULE

!~~~> use the memory buffers (tapes) implementation
  USE TLM_TAPES
!~~~> use the BLAS-like utility functions for vector computations
  USE UTILS
  IMPLICIT NONE
  PUBLIC
  SAVE

!~~~> maximum number of stages = s_Max = 16 (for DOPRI8(6) with dense output)
!~~~> A (s_Max by s_Max array, double precision),
!     b (s_Max by 1 array, double precision),
!     c (s_Max by 1 array, double precision) : coefficients for the Runge-Kutta methods
!~~~> e (s_Max by 1 array, double precision) = error coefficients
!~~~> Butcher tableau :
!           c | A
!           __|___
!             | b
  INTEGER, PARAMETER :: s_Max = 16
  DOUBLE PRECISION :: A(s_Max,s_Max) 
  DOUBLE PRECISION :: c(s_Max), b(s_Max), e(s_Max)
!~~~> dense output coefficients for the DOPRI8(6) method
  DOUBLE PRECISION :: d(s_Max,s_Max)
!~~~>  Statistics on the work performed by the Runge-Kutta method
!~~~> these are the indices of the ISTATUS and RSTATUS vectors returned by RKINT
  INTEGER, PARAMETER :: Nrhs=1, & !ISTATUS(Nrhs) = number of forward model ODE right-hand side evaluations
  			Nstp=2, & !ISTATUS(Nstp) = total number of time steps taken
  			Nacc=3, & !ISTATUS(Nacc) = total number of accepted time steps
                        Nrej=4, & !ISTATUS(Nrej) = total number of rejected time steps
                        Nrhs_tlm=5, & !ISTATUS(Nrhs_tlm) = number of TLM ODE right-hand side evaluations
                        Ntexit_tlm=1, & !RSTATUS(Ntexit_tlm) = Texit, the time corresponding to the
                                    !computed TLM solution upon return
                        Nhexit_tlm=2, & !RSTATUS(Nhexit_tlm) = Hexit, last accepted step before exit
                        Nhnew_tlm=3 !RSTATUS(Nhnew_tlm) = Hnew, last predicted step (not yet taken)

 CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE RKINT_TLM(N,Y,RHS,DY,TLM_RHS,               &
           Tstart,Tend,                                &
           Nd, Nc,                                     &
           AbsTol, RelTol,                             &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!    Solves the system y'=F(t,y) using an explicit Runge-Kutta method.
!    The user can choose from a range of options the desired RK method.
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
!     N  (integer)      = size of the ODE system
!     RHS (external)    = name (external) of subroutine computing the value of F(T,Y)
!                  SUBROUTINE RHS(N,T,Y,F)
!                  DOUBLE PRECISION T, Y(N),F(N)
!                  F(1) = .... etc.
!
!     TLM_RHS (external) = name (external) of the subroutine computing the right-hand-side
!                of the tangent linear ODE. It is recommended that the user
!                employ the forward mode of automatic differentiation to generate
!                the code for this subroutine (this allows a cheap evaluation of
!                the Jacobian-vector product (dF/dy)*DY as well as of the derivative
!                dF/dp_i).
!                  SUBROUTINE TLM_RHS(N,T1,Y1,DY1,TLM1)
!                       **Inputs:**
!                   INTEGER :: N !forward model/TLM state size
!                   DOUBLE PRECISION :: T1 !current time moment
!                   DOUBLE PRECISION :: Y1(N), DY1(N) !forward model state and
!                                                     !TLM state, respectively (at T1)
!                       **Output:**
!                   DOUBLE PRECISION :: TLM1(N)  !TLM ODE RHS evaluated at T1
!
!-    Y(N) (N by 1 double precision array)
!                  = vector of initial conditions for the forward model (at T=Tstart)
!
!-    DY(N) (N by 1 double precision array)
!                  = vector of initial conditions for the TLM model (at T=Tstart)
!
!     Tstart  (double precision)
!     Tend    (double precision)
!-   [Tstart,Tend] = time range of integration; RKINT requires Tstart <= Tend
!
!-    Nd (integer)  = number of steps between two consecutive checkpoints,
!               i.e. a checkpoint is written every Nd integration steps
!               Precondition: if checkpointing is enabled, RKINT requires that Nd > 0
!
!-    RelTol (2*N by 1 array, double precision)
!     AbsTol (2*N by 1 array, double precision)
!                          = user precribed tolerances for the forward
!                            model solution AbsTol(1:N), RelTol(1:N) and 
!                            forward sensitivities AbsTol(N+1:2*N), RelTol(N+1:2*N)
!
!     ICNTRL (20 by 1 array, integer)
!-    ICNTRL(1:20)  =  integer input parameters
!
!     RCNTRL (20 by 1 array, double precision)
!-    RCNTRL(1:20)  =  real input parameters
!
!     For more information on ICNTRL and RCNTRL, please see the INPUT PARAMETERS
!     section below.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!~~~>     INPUT PARAMETERS:
!
!    Note: For input parameters equal to zero the default values of the
!       corresponding variables are used.
!
!    ICNTRL(1) = 0: AbsTol, RelTol are N-dimensional vectors
!              = 1: AbsTol, RelTol are scalars
!
!    ICNTRL(2)  -> maximum number of integration steps
!        For ICNTRL(2)=0) the default value of 100000 is used
!
!    ICNTRL(3)  -> selection of a particular (explicit) embedded
!                 RK method with error control
!        = 0 :    DOPRI5(4) (default)
!        = 1 :    2nd order RK
!        = 2 :    Heun's method 3(2)
!        = 3 :    RK4(3)
!        = 4 :    DOPRI5(4)
!        = 5 :    RK6(5)9FM
!        = 6 :    DOPRI8(6)
!
!    ICNTRL(4)  -> Checkpointing flag
!        = 0 :    No checkpointing (trajectory is stored in memory)
!        = 1 :    Enable checkpointing (a checkpoint is written to disk every
!                 Nd steps)
!
!    ICNTRL(5)  -> Enable/Disable checkpointing and memory buffer mechanisms
!        = 0 :    Disabled (use for FWD+TLM model integration test runs; no data is stored)
!        = 1 :    Enabled (store information for adjoint sensitivity analysis)
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
!    T (double precision)   -> time for which the solution has been computed
!                     (after successful return T=Tend).
!
!    Y(N) (N by 1 array, double precision)
!                  -> Numerical solution of the forward model at T
!
!    DY(N) (N by 1 array, double precision)
!                  -> Numerical solution of the tangent linear model at T
!
!    IERR (integer)   -> Reports on successfulness upon return:
!                    = 1 for success
!                    < 0 for error (value equals error code)
!
!    ISTATUS(1)  -> Total number of forward model ODE right hand side evaluations
!    ISTATUS(2)  -> Total number of steps taken
!    ISTATUS(3)  -> Total number of accepted steps
!    ISTATUS(4)  -> Total number of rejected steps (except at very beginning)
!    ISTATUS(5)  -> Total number of TLM ODE right hand side evaluations
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
   INTEGER,       INTENT(IN)          :: N
   INTEGER,       INTENT(IN)          :: Nd
   INTEGER,       INTENT(OUT)         :: Nc
   DOUBLE PRECISION, INTENT(INOUT)    :: Y(N), DY(N)
   DOUBLE PRECISION, INTENT(INOUT)    :: Tstart,Tend
   DOUBLE PRECISION, INTENT(IN)       :: AbsTol(2*N),RelTol(2*N)
   INTEGER,       INTENT(IN)          :: ICNTRL(20)
   DOUBLE PRECISION, INTENT(IN)       :: RCNTRL(20)
   INTEGER,       INTENT(INOUT)       :: ISTATUS(20)
   DOUBLE PRECISION, INTENT(INOUT)    :: RSTATUS(20)
   INTEGER, INTENT(OUT)               :: IERR
   EXTERNAL RHS, TLM_RHS

!~~~>  Local variables
   DOUBLE PRECISION :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   DOUBLE PRECISION :: Hmin, Hmax, Hstart
   DOUBLE PRECISION :: Texit
   INTEGER       :: i, UplimTol, Max_no_steps
   LOGICAL       :: VectorTol
   LOGICAL       :: Checkpointing
   LOGICAL       :: AdjSensitivity

!~~~>  Double precision constants
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.d0, ONE  = 1.d0
!~~~>   Runge-Kutta method codes
   INTEGER, PARAMETER :: RK2=1, RK3=2, RK4=3, RK5=4, RK6=5, RK8=6
!~~~> this is used as the initial time step if no time step is provided in RCNTRL(3)
   DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0d-5

!~~~> rk_S = number of stages for the user-selected Runge-Kutta method
!~~~> rkMethod = code of the user-selected Runge-Kutta method
   INTEGER :: rk_S, rkMethod
!~~~> rkELO = order of accuracy for the user-selected RK method
! used in the time step and error control mechanism
   DOUBLE PRECISION :: rk_ELO

!~~~>  Initialize statistics
   ISTATUS(1:8) = 0
   RSTATUS(1:3) = ZERO

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!~~~>  For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:2*N) and RelTol(1:2*N)
   IF (ICNTRL(1) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = 2*N
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   All errors in the input parameters are handled by the rk_ErrorMsg subroutine

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(2) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(2) > 0) THEN
      Max_no_steps=ICNTRL(2)
   ELSE
      !error: the maximum number of steps cannot be negative!
      PRINT * ,'User-selected max no. of steps: ICNTRL(2)=',ICNTRL(2)
      CALL rk_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>   Initialize the particular Runge-Kutta method selected
   SELECT CASE (ICNTRL(3))
     CASE (1)
       CALL RK23
     CASE (2)
       CALL RK3_Heun
     CASE (3)
       CALL RK43
     CASE (0,4)
       CALL Dopri54
     CASE (5)
       CALL RK659FM
     CASE (6)
       CALL Dopri86
     CASE DEFAULT
       !invalid Runge-Kutta method selection
       PRINT * , 'Unknown Runge-Kutta method: ICNTRL(3)=',ICNTRL(3) 
       CALL rk_ErrorMsg(-9,Tstart,ZERO,IERR)
       RETURN
   END SELECT

   rkMethod = ICNTRL(3)

!~~~> Determine if file checkpointing is enabled
!~~~> If file checkpointing is enabled, RKINT requires that Nd > 0
   IF (ICNTRL(4) == 0) THEN
      !No checkpointing
      Checkpointing = .FALSE.
   ELSEIF (ICNTRl(4) == 1) THEN
      Checkpointing = .TRUE.
      IF (Nd .LE. 1) THEN
         !error: user enabled checkpointing but specified Nd <= 0
         PRINT *, 'Invalid value for parameter Nd = ', Nd
         CALL rk_ErrorMsg(-11,Tstart,ZERO,IERR)
      END IF
   ELSE
      !error: user specified an invalid checkpointing option
      PRINT * ,'Invalid checkpointing option: ICNTRL(4)=',ICNTRL(4)
      CALL rk_ErrorMsg(-10,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~> Determine if the user wants to enable the memory buffer mechanism
!~~~> Buffering must be enabled for adjoint sensitivity analysis studies
   IF (ICNTRL(5) == 0) THEN
      !Only a tangent linear model integration will be performed; no values will be stored
      !in memory and no file checkpoints will be written. This option takes precedence
      !over the file checkpointing option specified by ICNTRL(4)
      AdjSensitivity = .FALSE.
   ELSEIF (ICNTRL(5) == 1) THEN
      !Memory buffering is enabled
      AdjSensitivity = .TRUE.
   ELSE
      !Invalid option
      PRINT * ,'Invalid option: ICNTRL(5)=',ICNTRL(5)
      CALL rk_ErrorMsg(-12,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Unit roundoff or machine epsilon (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      !user did not impose a nonzero lower bound on the step size
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      !set lower bound on the step size as requested
      Hmin = RCNTRL(1)
   ELSE
      !error: lower bound on the step size cannot be negative
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL rk_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      !user did not impose an upper bound on the time step size
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      !limit Hmax as requested
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      !the upper limit for the time step size cannot be negative
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL rk_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      !User did not specify a starting step size
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      !User did specify a starting step size
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL rk_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      !set default value for FacMin
      FacMin = 0.2d0
   ELSEIF (RCNTRL(4) > ZERO) THEN
      !set user-requested value for FacMin
      FacMin = RCNTRL(4)
   ELSE
      !invalid value for FacMin: cannot be smaller than zero
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL rk_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(5) == ZERO) THEN
      !set default value for FacMax
      FacMax = 5.d0
   ELSEIF (RCNTRL(5) > ZERO) THEN
      !set user-requested value for FacMax
      FacMax = RCNTRL(5)
   ELSE
      !invalid value for FacRej: cannot be smaller than or equal to zero
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL rk_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      !set default value for FacRej
      FacRej = 0.1d0
   ELSEIF (RCNTRL(6) > ZERO) THEN
      !set user-requested value for FacRej
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL rk_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
   !default value for FacSafe
      FacSafe = 0.9d0
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      !invalid value for FacSafe: cannot be smaller than zero
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL rk_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

!~~~>  Check if tolerances are 'reasonable', i.e.:
!   0.d0 < AbsTol(i)
!   10.d0*Roundoff < RelTol(i) < 1.d0 for every i = 1..N
   DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.d0*Roundoff) &
         .OR. (RelTol(i) >= 1.0) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL rk_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
   END DO

!~~~>   Make sure that Tstart <= Tend - this is required by RKINT_TLM
   IF (Tstart .gt. Tend) THEN
       PRINT * ,'User-selected Tstart > Tend.'
       CALL rk_ErrorMsg(-8,Tstart,ZERO,IERR)
       RETURN
   END IF

!~~~>  CALL Runge-Kutta core integrator
   CALL RKTlm_Integrator(Y, RHS, DY,             &
        TLM_RHS,                                 &
        Nd, Nc,                                  &
        Tstart, Tend, Texit,                     &
        AbsTol, RelTol,                          &
!  Integration parameters
        VectorTol, Max_no_steps,                 &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
        Checkpointing,                           &
        AdjSensitivity,                          &
!  Real and integer output parameters
        RSTATUS, ISTATUS,                        &
!  Error indicator
        IERR)

 CONTAINS


SUBROUTINE RK23
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Initializes coefficients (the Butcher tableau) for the RK2(3) 
! Runge-Kutta method. This is a 2(3) Fehlberg-type method: second
! order accurate with a third order error control mechanism
! 
! A = [ 0 0 0; 
!       1 0 0; 
!       0.25 0.25 0 ]
!
! b = [ 0.5 0.5 0 ]
!
! c = [ 0 1 0.5 ]
!
! e = [ 1/3 1/3 -2/3 ] 
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  IMPLICIT NONE

  rk_S = 3
  rk_ELO = 2.d0
  rkMethod = RK2

  c = 0.d0
  b = 0.d0
  A = 0.d0

  c(2) = 1.d0
  c(3) = 0.5d0

  A(2,1) = 1.d0
  A(3,1) = 0.25d0
  A(3,2) = 0.25d0

  b(1) = 0.5d0
  b(2) = 0.5d0

  e = 0.d0
  e(1) = 1.d0/3.d0
  e(2) = 1.d0/3.d0
  e(3) = -2.d0/3.d0

END SUBROUTINE RK23



SUBROUTINE RK3_Heun
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Initializes parameters (Butcher tableau) for the RK3(2)
! Runge-Kutta method built on top of Heun's third order method.
! It incorporates a second-order accurate error control mechanism.
!
! A = [ 0 0 0 0;
!       1/3 0 0 0; 
!       0 2/3 0 0;
!       0.25 0 0.75 0]
!
! b = [ 0.25 0 0.75 0 ]
!
! c = [ 0 1/3 2/3 0 ]
!
! e = [ 0.25 -0.5 0.25 0 ]
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  rk_S = 4
  rk_ELO = 3.d0
  rkMethod = RK3

  c = 0
  c(2) = 1.d0 / 3.d0
  c(3) = 2.d0 / 3.d0

  A = 0.d0
  A(2,1) = 1.d0 / 3.d0
  A(3,2) = 2.d0 / 3.d0
  A(4,1) = 0.25d0
  A(4,3) = 0.75d0

  b = 0.d0
  b(2) = 0.5d0
  b(3) = 0.5d0

  e  = A(4,:) - b
  b = A(4,:)

END SUBROUTINE RK3_Heun



SUBROUTINE RK43
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Initializes coefficients (the Butcher tableau) for the RK4(3)
! Runge-Kutta method. This method is built by appending an
! additional stage to the well-known "3/8-rule", thus incorporating
! a 3rd order error control mechanism.
!
! A = [ 0 0 0 0 0;
!       1/3 0 0 0 0;
!       -1/3 1 0 0 0;
!       1 -1 -1 0 0;
!       1/8 3/8 3/8 1/8 0];
!
! b = [ 1/8 3/8 3/8 1/8 0]
! 
! c = [ 0 1/3 2/3 1 1 ]
!
! e = [ 1/24 -1/8 1/8 1/8 -1/6 ]
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  rk_S = 5
  rk_ELO = 4.d0
  rkMethod = RK4

  b = 0.d0
  c = 0.d0
  A = 0.d0

  c(2) = 1.d0 / 3.d0
  c(3) = 2.d0 / 3.d0
  c(4) = 1.d0
  c(5) = 1.d0

  A(2,1) = 1.d0 / 3.d0
  A(3,1) = -1.d0 / 3.d0
  A(3,2) = 1.d0
  A(4,1) = 1.d0
  A(4,2) = -1.d0
  A(4,3) = 1.d0
  A(5,1) = 0.125d0
  A(5,2) = 0.375d0
  A(5,3) = 0.375d0
  A(5,4) = 0.125d0

  b(1) = 1.d0 / 12.d0
  b(2) = 0.5d0
  b(3) = 0.25d0
  b(5) = 1.d0 / 6.d0

  e = A(5,:) - b
  b = A(5,:)

END SUBROUTINE RK43



SUBROUTINE Dopri54
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Initializes parameters (Butcher tableau) for the DOPRI5(4) 
! Runge-Kutta method. For further details on the method coeffcients
! please see "A family of embedded Runge-Kutta formulae" by J. Dormand
! and P.J. Prince (J. Comput. Appl. Math., vol.6, 1980, p. 19-26)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  IMPLICIT NONE

  !number of stages
  rk_S = 7
  rk_ELO = 5.d0
  rkMethod = RK5

  c = 0.d0
  !set up the coefficients
  c(2) = 0.2d0
  c(3) = 0.3d0
  c(4) = 0.8d0
  c(5) = 8.d0/9.d0
  c(6) = 1.d0
  c(7) = 1.d0

  A = 0.d0
  A(2,1) = 0.2d0
  A(3,1) = 3.d0 / 40.d0
  A(3,2) = 9.d0 / 40.d0
  A(4,1) = 44.d0 / 45.d0
  A(4,2) = -56.d0 / 15.d0
  A(4,3) = 32.d0 / 9.d0
  A(5,1) = 19372.d0 / 6561.d0
  A(5,2) = -25360.d0 / 2187.d0
  A(5,3) = 64448.d0 / 6561.d0
  A(5,4) = -212.d0 / 729.d0
  A(6,1) = 9017.d0 / 3168.d0
  A(6,2) = -355.d0 / 33.d0
  A(6,3) = 46732.d0 / 5247.d0
  A(6,4) = 49.d0 / 176.d0
  A(6,5) = -5103.d0 / 18656.d0
  A(7,1) = 35.d0 / 384.d0
  A(7,3) = 500.d0 / 1113.d0
  A(7,4) = 125.d0 / 192.d0
  A(7,5) = -2187.d0 / 6784.d0
  A(7,6) = 11.d0 / 84.d0

  b = A(7,:)

  e = 0.d0
  e(1) = 71.d0/57600.D0
  e(3) = -71.d0/16695.d0
  e(4) = 71.d0/1920.d0
  e(5) = -17253.d0/339200.d0
  e(6) = 22.d0/525.d0
  e(7) = -1.d0/40.d0

END SUBROUTINE Dopri54



SUBROUTINE RK659FM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Initializes parameters (Butcher tableau) for the Runge-Kutta 6(5) 
! method with 5th order dense output built on top of RK6(5)9FM (see
! "Global error estimation with Runge-Kutta triples" by 
! J. R. Dormand and M. A. Lockyer and N. E. McGorrigan and P. J. Prince,
! Comput. Math. Appl., 1989, no. 9, vol. 18, p. 835-846)
! with 6th order dense output (see "Continuous approximation with embedded 
! Runge-Kutta methods", T. S. Baker and J. R. Dormand and J. P. Gilmore 
! and P. J. Prince, Appl. Numer. Math., 1996, vol 22, p. 51-62)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  IMPLICIT NONE

  !number of stages
  rk_S = 12
  rk_ELO = 6.d0
  rkMethod = RK6

  c = 0.d0
  c(2) = 1.d0/9.d0
  c(3) = 1.d0/6.d0
  c(4) = 0.25d0
  c(5) = 5.d0/9.d0
  c(6) = 0.5d0
  c(7) = 48.d0/49.d0
  c(8) = 1.d0
  c(9) = 1.d0
  c(10) = 0.5d0
  c(11) = 1.d0/6.d0
  c(12) = 5.d0/12.d0

  A = 0.d0
  A(2,1) = 1.d0/9.d0
  A(3,1) = 1.d0/24.d0
  A(3,2) = 1.d0/8.d0
  A(4,1) = 1.d0/16.d0
  A(4,3) = 3.d0/16.d0
  A(5,1) = 280.d0/729.d0
  A(5,3) = -325.d0/243.d0
  A(5,4) = 1100.d0/729.d0
  A(6,1) = 6127.d0/14680.d0
  A(6,3) = -1077.d0/734.d0
  A(6,4) = 6494.d0/4037.d0
  A(6,5) = -9477.d0/161480.d0
  A(7,1) = -13426273320.d0/14809773769.d0
  A(7,3) = 4192558704.d0/2115681967.d0
  A(7,4) = 14334750144.d0/14809773769.d0
  A(7,5) = 117092732328.d0/14809773769.d0
  A(7,6) = -361966176.d0/40353607.d0
  A(8,1) = -2340689.d0/1901060.d0
  A(8,3) = 31647.d0/13579.d0
  A(8,4) = 253549596.d0/149518369.d0
  A(8,5) = 10559024082.d0/977620105.d0
  A(8,6) = -152952.d0/12173.d0
  A(8,7) = -5764801.d0/186010396.d0
  A(9,1) = 203.d0/2880.d0
  A(9,4) = 30208.d0/70785.d0
  A(9,5) = 177147.d0/164560.d0
  A(9,6) = -536.d0/705.d0
  A(9,7) = 1977326743.d0/3619661760.d0
  A(9,8) = -259.d0/720.d0
  A(10,1) = 33797.d0/460800.d0
  A(10,4) = 27757.d0/70785.d0
  A(10,5) = 7923501.d0/26329600.d0
  A(10,6) = -927.d0/3760.d0
  A(10,7) = -3314760575.d0/23165835264.d0
  A(10,8) = 2479.d0/23040.d0
  A(10,9) = 1.d0/64.d0
  A(11,1) = 5843.d0/76800.d0
  A(11,4) = 464.d0/2673.d0
  A(11,5) = 353997.d0/1196800.d0
  A(11,6) = -15068.d0/57105.d0
  A(11,7) = -282475249.d0/3644974080.d0
  A(11,8) = 8678831.d0/156245760.d0
  A(11,9) = 116113.d0/11718432.d0
  A(11,10) = -25.d0/243.d0
  A(12,1) = 15088049.d0/199065600.d0
  A(12,4) = 0.4d0
  A(12,5) = 92222037.d0/268083200.d0
  A(12,6) = -433420501.d0/1528586640.d0
  A(12,7) = -11549242677007.d0/83630285291520.d0
  A(12,8) = 2725085557.d0/26167173120.d0
  A(12,9) = 235429367.d0/16354483200.d0
  A(12,10) = -90924917.d0/1040739840.d0
  A(12,11) = -271149.d0/21414400.d0

  b = 0.d0
  b(1) = 36567.d0/458800.d0
  b(4) = 9925984.d0/27063465.d0
  b(5) = 85382667.d0/117968950.d0
  b(6) = -310378.d0/808635.d0
  b(7) = 262119736669.d0/345979336560.d0
  b(8) = -0.5d0
  b(9) = -101.d0/2294.d0

  e = A(9,:) - b

END SUBROUTINE RK659FM



SUBROUTINE Dopri86
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Initializes parameters (Butcher tableau and dense output
! coefficients) for the DOPRI8(6) Runge-Kutta method with 7th 
! order dense output. For further information on the method 
! and dense output coefficients, please see 
! "High order embedded Runge-Kutta formulae" by P.J.Prince and 
! J. Dormand, J. Comput Appl Math, 1981, vol. 7, p. 67-76, as well
! as the DOPRI853 subroutine by E. Hairer and G. Wanner (source code
! to be found in the appendix of "Solving Ordinary Differential 
! Equations: Nonstiff Problems", Springer-Verlag, 1993)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  IMPLICIT NONE

  rk_S = 16
  rk_ELO = 8.d0
  rkMethod = RK8

  c = 0.d0
  c(2)  = 0.526001519587677318785587544488D-01
  c(3)  = 0.789002279381515978178381316732D-01
  c(4)  = 0.118350341907227396726757197510D+00
  c(5)  = 0.281649658092772603273242802490D+00
  c(6)  = 0.333333333333333333333333333333D+00
  c(7)  = 0.25D+00
  c(8)  = 0.307692307692307692307692307692D+00
  c(9)  = 0.651282051282051282051282051282D+00
  c(10) = 0.6D+00
  c(11) = 0.857142857142857142857142857142D+00
  c(12) = 1.D0
  c(13) = 1.D0
  c(14) = 0.1D+00
  c(15) = 0.2D+00
  c(16) = 0.777777777777777777777777777778D+00

  b = 0.D0
  b(1) =   5.42937341165687622380535766363D-2
  b(6) =   4.45031289275240888144113950566D0
  b(7) =   1.89151789931450038304281599044D0
  b(8) =  -5.8012039600105847814672114227D0
  b(9) =   3.1116436695781989440891606237D-1
  b(10) = -1.52160949662516078556178806805D-1
  b(11) =  2.01365400804030348374776537501D-1
  b(12) =  4.47106157277725905176885569043D-2

  A = 0.D0
  A(2,1) =    5.26001519587677318785587544488D-2
  A(3,1) =    1.97250569845378994544595329183D-2
  A(3,2) =    5.91751709536136983633785987549D-2
  A(4,1) =    2.95875854768068491816892993775D-2
  A(4,3) =    8.87627564304205475450678981324D-2
  A(5,1) =    2.41365134159266685502369798665D-1
  A(5,3) =   -8.84549479328286085344864962717D-1
  A(5,4) =    9.24834003261792003115737966543D-1
  A(6,1) =    3.7037037037037037037037037037D-2
  A(6,4) =    1.70828608729473871279604482173D-1
  A(6,5) =    1.25467687566822425016691814123D-1
  A(7,1) =    3.7109375D-2
  A(7,4) =    1.70252211019544039314978060272D-1
  A(7,5) =    6.02165389804559606850219397283D-2
  A(7,6) =   -1.7578125D-2

  A(8,1) =    3.70920001185047927108779319836D-2
  A(8,4) =    1.70383925712239993810214054705D-1
  A(8,5) =    1.07262030446373284651809199168D-1
  A(8,6) =   -1.53194377486244017527936158236D-2
  A(8,7) =    8.27378916381402288758473766002D-3
  A(9,1) =    6.24110958716075717114429577812D-1
  A(9,4) =   -3.36089262944694129406857109825D0
  A(9,5) =   -8.68219346841726006818189891453D-1
  A(9,6) =    2.75920996994467083049415600797D1
  A(9,7) =    2.01540675504778934086186788979D1
  A(9,8) =   -4.34898841810699588477366255144D1
  A(10,1) =   4.77662536438264365890433908527D-1
  A(10,4) =  -2.48811461997166764192642586468D0
  A(10,5) =  -5.90290826836842996371446475743D-1
  A(10,6) =   2.12300514481811942347288949897D1
  A(10,7) =   1.52792336328824235832596922938D1
  A(10,8) =  -3.32882109689848629194453265587D1
  A(10,9) =  -2.03312017085086261358222928593D-2

  A(11,1) =  -9.3714243008598732571704021658D-1
  A(11,4) =   5.18637242884406370830023853209D0
  A(11,5) =   1.09143734899672957818500254654D0
  A(11,6) =  -8.14978701074692612513997267357D0
  A(11,7) =  -1.85200656599969598641566180701D1
  A(11,8) =   2.27394870993505042818970056734D1
  A(11,9) =   2.49360555267965238987089396762D0
  A(11,10) = -3.0467644718982195003823669022D0
  A(12,1) =   2.27331014751653820792359768449D0
  A(12,4) =  -1.05344954667372501984066689879D1
  A(12,5) =  -2.00087205822486249909675718444D0
  A(12,6) =  -1.79589318631187989172765950534D1
  A(12,7) =   2.79488845294199600508499808837D1
  A(12,8) =  -2.85899827713502369474065508674D0
  A(12,9) =  -8.87285693353062954433549289258D0
  A(12,10) =  1.23605671757943030647266201528D1
  A(12,11) =  6.43392746015763530355970484046D-1

  A(13,:) = b

  A(14,1) =  5.61675022830479523392909219681D-2
  A(14,7) =  2.53500210216624811088794765333D-1
  A(14,8) = -2.46239037470802489917441475441D-1
  A(14,9) = -1.24191423263816360469010140626D-1
  A(14,10) =  1.5329179827876569731206322685D-1
  A(14,11) =  8.20105229563468988491666602057D-3
  A(14,12) =  7.56789766054569976138603589584D-3
  A(14,13) = -8.298D-3

  A(15,1) =  3.18346481635021405060768473261D-2
  A(15,6) =  2.83009096723667755288322961402D-2
  A(15,7) =  5.35419883074385676223797384372D-2
  A(15,8) = -5.49237485713909884646569340306D-2
  A(15,11) = -1.08347328697249322858509316994D-4
  A(15,12) =  3.82571090835658412954920192323D-4
  A(15,13) = -3.40465008687404560802977114492D-4
  A(15,14) =  1.41312443674632500278074618366D-1
  A(16,1) = -4.28896301583791923408573538692D-1
  A(16,6) = -4.69762141536116384314449447206D0
  A(16,7) =  7.68342119606259904184240953878D0
  A(16,8) =  4.06898981839711007970213554331D0
  A(16,9) =  3.56727187455281109270669543021D-1
  A(16,13) = -1.39902416515901462129418009734D-3
  A(16,14) =  2.9475147891527723389556272149D0
  A(16,15) = -9.15095847217987001081870187138D0

  e = 0.d0
  e(1) =  0.1312004499419488073250102996D-01
  e(6) = -0.1225156446376204440720569753D+01
  e(7) = -0.4957589496572501915214079952D+00
  e(8) =  0.1664377182454986536961530415D+01
  e(9) = -0.3503288487499736816886487290D+00
  e(10) =  0.3341791187130174790297318841D+00
  e(11) =  0.8192320648511571246570742613D-01
  e(12) = -0.2235530786388629525884427845D-01

!~~~> Coefficients used for the 7th order dense output
!~~~> For the code that inspired this implementation, see E. Hairer and G. Wanner's
!~~~> DOPRI853 routine at http://www.unige.ch/~hairer/prog/nonstiff/dop853.f
  d = 0.d0
  d(4,1)  = -0.84289382761090128651353491142D+01
  d(4,6)  =  0.56671495351937776962531783590D+00
  d(4,7)  = -0.30689499459498916912797304727D+01
  d(4,8)  =  0.23846676565120698287728149680D+01
  d(4,9)  =  0.21170345824450282767155149946D+01
  d(4,10) = -0.87139158377797299206789907490D+00
  d(4,11) =  0.22404374302607882758541771650D+01
  d(4,12) =  0.63157877876946881815570249290D+00
  d(4,13) = -0.88990336451333310820698117400D-01
  d(4,14) =  0.18148505520854727256656404962D+02
  d(4,15) = -0.91946323924783554000451984436D+01
  d(4,16) = -0.44360363875948939664310572000D+01

  d(5,1)  =  0.10427508642579134603413151009D+02
  d(5,6)  =  0.24228349177525818288430175319D+03
  d(5,7)  =  0.16520045171727028198505394887D+03
  d(5,8)  = -0.37454675472269020279518312152D+03
  d(5,9)  = -0.22113666853125306036270938578D+02
  d(5,10) =  0.77334326684722638389603898808D+01
  d(5,11) = -0.30674084731089398182061213626D+02
  d(5,12) = -0.93321305264302278729567221706D+01
  d(5,13) =  0.15697238121770843886131091075D+02
  d(5,14) = -0.31139403219565177677282850411D+02
  d(5,15) = -0.93529243588444783865713862664D+01
  d(5,16) =  0.35816841486394083752465898540D+02

  d(6,1) =  0.19985053242002433820987653617D+02
  d(6,6) = -0.38703730874935176555105901742D+03
  d(6,7) = -0.18917813819516756882830838328D+03
  d(6,8) =  0.52780815920542364900561016686D+03
  d(6,9) = -0.11573902539959630126141871134D+02
  d(6,10) =  0.68812326946963000169666922661D+01
  d(6,11) = -0.10006050966910838403183860980D+01
  d(6,12) =  0.77771377980534432092869265740D+00
  d(6,13) = -0.27782057523535084065932004339D+01
  d(6,14) = -0.60196695231264120758267380846D+02
  d(6,15) =  0.84320405506677161018159903784D+02
  d(6,16) =  0.11992291136182789328035130030D+02

  d(7,1)  = -0.25693933462703749003312586129D+02
  d(7,6)  = -0.15418974869023643374053993627D+03
  d(7,7)  = -0.23152937917604549567536039109D+03
  d(7,8)  =  0.35763911791061412378285349910D+03
  d(7,9)  =  0.93405324183624310003907691704D+02
  d(7,10) = -0.37458323136451633156875139351D+02
  d(7,11) =  0.10409964950896230045147246184D+03
  d(7,12) =  0.29840293426660503123344363579D+02
  d(7,13) = -0.43533456590011143754432175058D+02
  d(7,14) =  0.96324553959188282948394950600D+02
  d(7,15) = -0.39177261675615439165231486172D+02
  d(7,16) = -0.14972683625798562581422125276D+03

END SUBROUTINE Dopri86



SUBROUTINE rk_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages. 
! Inputs:
!      Code (integer) : Error code
!      T (double precision) : Time at which the error occured
!      H (double precision) : Last time step taken before the error occured
! Output:
!      IERR (integer) : Error indicator. Upon exit from rk_ErrorMsg, IERR = Code
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR

   IERR = Code
   PRINT * , &
     'Forced exit from RKINT_TLM due to the following error:' 

   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)
      PRINT * , '--> Tapes are not yet allocated. Could not store forward integration data.'
    CASE (-3)
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4) 
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)
      PRINT * , '--> User selected Tstart > Tend. ', &
            'Backward time integration in the forward mode is not supported.'
    CASE (-9)
      PRINT * , '--> Selected Runge-Kutta method not implemented'
    CASE (-10)
      PRINT * , '--> Improper checkpointing option.'
    CASE (-11)
      PRINT * , '--> Improper value for Nd. Nd must be greater than zero when checkpointing is enabled.'
    CASE (-12)
      PRINT * , '--> Improper option for memory buffers + checkpointing mechanisms. ', & 
            'Must choose enabled (1) or disabled (0).'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT

   PRINT *, "T=", T, "and H=", H

END SUBROUTINE rk_ErrorMsg



SUBROUTINE RKTLM_Integrator (Y, RHS, DY,           &
        TLM_RHS,                                   &
        Nd, Nc,                                    &
        Tstart, Tend, T,                           &
        AbsTol, RelTol,                            &
!~~~> Integration parameters
        VectorTol, Max_no_steps,                   &
        Roundoff, Hmin, Hmax, Hstart,              &
        FacMin, FacMax, FacRej, FacSafe,           &
        Checkpointing,                             &
        AdjSensitivity,                            &
        RSTATUS, ISTATUS,                          &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Core Runge-Kutta integrator. Performs the actual integration of the ODE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Input parameters:
!         Y (N by 1 array, double precision) : forward model state at Tstart
!         RHS (external) : pointer to (user-supplied) subroutine that evaluates
!                          the right-hand side of the forward model ODE
!         DY (N by 1 array, double precision) : TLM model state at Tstart
!         TLM_RHS (external) : pointer to (user-supplied) subroutine that evaluates
!                          the right-hand side of the TLM ODE
!         Nd (integer) : number of time steps between two consecutive checkpoints
!                        if Checkpointing == FALSE, the value of Nd is irrelevant
!         Tstart (double precision) : integration start time
!         Tend (double precision) : integration end time
!         AbsTol (N by 1 array, double precision) : absolute tolerances
!         RelTol (N by 1 array, double precision) : relative tolerances
!         VectorTol (logical) : vector tolerances (TRUE) or scalar tolerances (FALSE)
!         Max_no_steps (integer) : maximum number of integration steps
!         Roundoff (double precision) : machine epsilon
!         Hmin (double precision) : lower bound for the integration step size
!         Hmax (double precision) : upper bound for the integration step size
!         Hstart (double precision) : starting value for the integration step size
!         FacMin (double precision) : lower bound on step decrease factor
!         FacMax (double precision) : upper bound on step increase factor
!         FacRej (double precision) : step decrease factor after multiple rejections
!         FacSafe (double precision) : factor by which the new step is slightly smaller
!                                      than the predicted value
!         Checkpointing (logical) : enable (TRUE) / disable (FALSE) file checkpointing mechanism
!         AdjSensitivity (logical) : enable (TRUE) / disable (TRUE) memory buffering for
!                                    forward trajectory storage (needed for adjoint computations)
!
!
! Output parameters:
!
!         Nc (integer) : total number of file checkpoints written during the integration
!         RSTATUS (20 by 1 array, double precision) : integrator output parameters
!         ISTATUS (20 by 1 array, double precision) : integrator output parameters (statistics)
!
!         For more information on RSTATUS and ISTATUS, see the comments describing the RKINT subroutine
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   IMPLICIT NONE
   EXTERNAL RHS, TLM_RHS
!~~~> Input: the initial condition at Tstart; Output: the solution at T
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N)
!~~~> Input: the initial sensitivities at Tstart; Output: the TLM solution at T
   DOUBLE PRECISION, INTENT(INOUT) :: DY(N)
!~~~> Input: integration interval
   DOUBLE PRECISION, INTENT(INOUT) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   DOUBLE PRECISION, INTENT(OUT) ::  T
!~~~> Input: if checkpointing is enabled, the number of integration steps
!            between two consecutive checkpoints
   INTEGER, INTENT(IN) :: Nd
!~~~> Output: the number of checkpoints written during the integration
!             Nc = 0 if checkpointing is disabled
   INTEGER, INTENT(OUT) :: Nc
!~~~> Input: tolerances
   DOUBLE PRECISION, INTENT(IN) ::  AbsTol(2*N), RelTol(2*N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: VectorTol
   DOUBLE PRECISION, INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   LOGICAL, INTENT(IN) :: Checkpointing
   LOGICAL, INTENT(IN) :: AdjSensitivity
   DOUBLE PRECISION, INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   DOUBLE PRECISION, INTENT(INOUT) :: RSTATUS(20)
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables for the new solution and RHS values at the 1st and later stages
   DOUBLE PRECISION :: Ynew(N), Fcn0(N), Fcn(N)
   DOUBLE PRECISION :: DYnew(N), Tlm0(N), Tlm(N)
!K and Ktlm = the stages for the FWD and TLM models
   DOUBLE PRECISION :: K(N,rk_S)
   DOUBLE PRECISION :: Ktlm(N,rk_S)
!K_buf, Ktlm_buf = buffers for tape writes
   DOUBLE PRECISION :: K_buf(N,s_Max), Ktlm_buf(N,s_Max)

   DOUBLE PRECISION :: H, Hnew, Fac, Tau
   DOUBLE PRECISION :: Err, Yerr(N), DYerr(N)
   DOUBLE PRECISION :: FullSol(2*N), FullSolNew(2*N), FullErr(2*N)

   INTEGER :: j, istage
   INTEGER:: allocErr = 0 !flag will be set if an error
                          !is encountered in the tapes module
   LOGICAL :: RejectLastH, RejectMoreH
!~~~>  Local parameters
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.d0, ONE  = 1.d0
   DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0d-5
!~~~> Integration steps counter
   INTEGER :: chkStpCount

!~~~> zero out the parts of Kbuf and Ktlm_buf that will not be used
   DO j=rk_S+1,s_Max
       CALL SET2ZERO(N,K_buf(:,j))
       CALL SET2ZERO(N,Ktlm_buf(:,j))
   END DO

   IF (AdjSensitivity) THEN
       !Check to see if the buffers have been allocated
       CALL tlm_tapes_CheckAllocated(allocErr)
       !Buffers have not yet been allocated!
       !Signal this error
       IF (allocErr .NE. 1) THEN
           !tapes have not yet been allocated!
           !signal this error and return from subroutine
           CALL rk_ErrorMsg(-2,Tstart,ZERO,IERR)
           RETURN
       END IF
   END IF

!~~~>  Initial preparations
   Nc = 0
   chkStpCount = 0
   T = Tstart
   RSTATUS(Nhexit_tlm) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.d0*Roundoff) H = DeltaMin

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.

!~~~> Time loop begins below

TimeLoop: DO WHILE ((T-Tend)+Roundoff <= ZERO)

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL rk_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1d0*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL rk_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL RHS(N,T,Y,Fcn0)
   ISTATUS(Nrhs) = ISTATUS(Nrhs) + 1
   CALL TLM_RHS(N,T,Y,DY,Tlm0)
   ISTATUS(Nrhs_tlm) = ISTATUS(Nrhs_tlm) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO

!~~~>   Compute the stages
Stage: DO istage = 1, rk_S

      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         !Fcn = Fcn0
         CALL WCOPY(N,Fcn0,1,Fcn,1)
         !Tlm = Tlm0
         CALL WCOPY(N,Tlm0,1,Tlm,1)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSE

         !Ynew = Y
         CALL WCOPY(N,Y,1,Ynew,1)
         DO j = 1, istage-1
           !Ynew = Ynew + H*A(istage,j)*K(:,j)
           CALL WAXPY(N,H*A(istage,j),K(:,j),1,Ynew,1)
         END DO

         !DYnew = DY
         CALL WCOPY(N,DY,1,DYnew,1)
         DO j = 1, istage-1
           CALL WAXPY(N,H*A(istage,j),Ktlm(:,j),1,DYnew,1)
         END DO

         !evaluate RHS at time t=Tau and y=Ynew
         Tau = T + c(istage)*H
         CALL RHS(N,Tau,Ynew,Fcn)
         ISTATUS(Nrhs) = ISTATUS(Nrhs) + 1

         CALL TLM_RHS(N,Tau,Ynew,DYnew,Tlm)
         ISTATUS(Nrhs_tlm) = ISTATUS(Nrhs_tlm) + 1

       END IF

       !K(:,istage) = Fcn
       CALL WCOPY(N,Fcn,1,K(:,istage),1)
       !Ktlm(:,istage) = Tlm
       CALL WCOPY(N,Tlm,1,Ktlm(:,istage),1)

   END DO Stage

!~~~>  Compute the new solution
   !Ynew = Y
   CALL WCOPY(N,Y,1,Ynew,1)
   DO j=1,rk_S
        !Ynew = Ynew + H*b(j)*K(:,j)
        CALL WAXPY(N,H*b(j),K(:,j),1,Ynew,1)
   END DO

   !DYnew = DY
   CALL WCOPY(N,DY,1,DYnew,1)
   DO j=1,rk_S
        !DYnew = DYnew + H*b(j)*Ktlm(:,j)
        CALL WAXPY(N,H*b(j),Ktlm(:,j),1,DYnew,1)
   END DO

!~~~>  Compute the error estimation
   !Yerr = 0.d0
   CALL SET2ZERO(N,Yerr)
   DO j=1,rk_S
        !Yerr = Yerr + H*e(j)*K(:,j)
        CALL WAXPY(N,H*e(j),K(:,j),1,Yerr,1)
   END DO

   !DYerr = 0.d0
   CALL SET2ZERO(N,DYerr)
   DO j=1,rk_S
        !DYerr = DYerr + H*e(j)*Ktlm(:,j)
        CALL WAXPY(N,H*e(j),Ktlm(:,j),1,DYerr,1)
   END DO

!~~~> Copy the FWD and TLM current state, new state and error vectors
!~~~> to expanded vectors (of size 2*N) used to compute the scaled error norm
   CALL WCOPY(N,Y,1,FullSol(1:N),1)
   CALL WCOPY(N,DY,1,FullSol(N+1:2*N),1)

   CALL WCOPY(N,Ynew,1,FullSolNew(1:N),1)
   CALL WCOPY(N,DYnew,1,FullSolNew(N+1:2*N),1)

   CALL WCOPY(N,Yerr,1,FullErr(1:N),1)
   CALL WCOPY(N,DYerr,1,FullErr(N+1:2*N),1)

!~~~> Compute the scaled norm of the error
   Err = rk_ErrorNorm ( 2*N, FullSol, FullSolNew, Fullerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/rk_ELO)))
   Hnew = H*Fac

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1

   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step

      IF (AdjSensitivity) THEN

        IF (Checkpointing) THEN

          !~~~> Checkpointing is enabled; check if it is time to write a file checkpoint
          IF (mod(ISTATUS(Nacc),Nd) .EQ. 0) THEN

            !~~~> Write checkpoint; store enough data to allow a hot restart
            CALL tlm_tapes_WriteCheckpt(T,H,Y,DY)
            !~~~> Increment number of checkpoints & flush tapes
            CALL tlm_tapes_Flush
            !~~~> Increment checkpoint counter
            Nc = Nc + 1

          END IF

          chkStpCount = chkStpCount + 1
        END IF

        !~~~>  Record the timestep, solution, current time and stages
        !K_buf(:,1:rk_S) = K
        !Ktlm_buf(:,1:rk_S) = Ktlm
        DO j=1,rk_S
            CALL WCOPY(N,K(:,j),1,K_buf(:,j),1)
            CALL WCOPY(N,Ktlm(:,j),1,Ktlm_buf(:,j),1)
        END DO
        CALL tlm_tapes_Write(T,H,K_buf,Ktlm_buf,Fcn0,Tlm0,Y,DY)

      END IF

      !increment the "accepted steps" counter
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      !update the solution
      !Y = Ynew
      CALL WCOPY(N,Ynew,1,Y,1)
      !DY = DYnew
      CALL WCOPY(N,DYnew,1,DY,1)
      !update current time
      T = T + H
      !make sure Hmin <= Hnew <= Hmax
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H)
      END IF
      RSTATUS(Nhexit_tlm) = H
      RSTATUS(Nhnew_tlm)  = Hnew
      RSTATUS(Ntexit_tlm) = T
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      !final update update step size
      H = Hnew
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED

   ELSE           !~~~> Reject step

      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1

   END IF ! Err <= 1

   END DO UntilAccepted

   END DO TimeLoop

   IF (AdjSensitivity) THEN

      CALL RHS(N,Tend,Y,Fcn0)
      ISTATUS(Nrhs) = ISTATUS(Nrhs) + 1

      CALL TLM_RHS(N,Tend,Y,DY,Tlm0)
      ISTATUS(Nrhs_tlm) = ISTATUS(Nrhs_tlm) + 1

      !save last FWD&TLM state (this is needed in the continuous adjoint run)
      !K_buf(:,1:rk_S) = K
      !Ktlm_buf(:,1:rk_S) = Ktlm
      DO j=1,rk_S
           CALL WCOPY(N,K(:,j),1,K_buf(:,j),1)
           CALL WCOPY(N,Ktlm(:,j),1,Ktlm_buf(:,j),1)
      END DO
      CALL tlm_tapes_Write(Tend,ZERO,K,Ktlm,Fcn0,Tlm0,Y,DY)

   END IF

!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

END SUBROUTINE RKTLM_Integrator



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


!End of TLM integrator routine
END SUBROUTINE RKINT_TLM



SUBROUTINE rktlm_AllocateTapes(CheckptEnabled, MaxSize, N)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Tape allocator -> calls tlm_tapes_Allocate
!~~~> This is the routine that the user should call from the driver program
!to allocate memory for the buffer mechanism
! INPUTS:
!       LOGICAL CheckptEnabled --> Checkpoint mechanism enabled/disabled
!                     = FALSE : Disabled (memory storage only)
!                     = TRUE  : Enabled (disk checkpoints + memory buffers)
!       INTEGER MaxSize --> Memory buffer capacity. This should be greater than
!                           or equal to the number of time steps between two
!                           consecutive checkpoints (if CheckptEnabled == TRUE) or
!                           the maximum number of integration time steps (if
!                           CheckptEnabled == FALSE).
!       INTEGER N --> Size of the TLM state
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: N
  LOGICAL, INTENT(INOUT) :: CheckptEnabled
  INTEGER, INTENT(INOUT) :: MaxSize

  CALL tlm_tapes_Allocate(CheckptEnabled, MaxSize, N, s_Max)

END SUBROUTINE rktlm_AllocateTapes



SUBROUTINE rktlm_DeallocateTapes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Tape deallocator -> calls tlm_tapes_Deallocate
!~~~> This subroutine is callable from the driver program
!~~~> Should be called after using the stored forward model trajectory to free
! all allocated memory buffers.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  IMPLICIT NONE

  CALL tlm_tapes_Deallocate

END SUBROUTINE rktlm_DeallocateTapes


!END OF MODULE
END MODULE RKTLM_MODULE