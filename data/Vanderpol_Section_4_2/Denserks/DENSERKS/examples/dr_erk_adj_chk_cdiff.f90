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


MODULE PERTURBATION

  IMPLICIT NONE
  PUBLIC
  SAVE

  DOUBLE PRECISION, PARAMETER :: dp1 = 0.001d0
  DOUBLE PRECISION, PARAMETER :: dp2 = 0.0005d0

END MODULE PERTURBATION

! Convection-diffusion problem parameters
MODULE CDIFF_PARAMS

  IMPLICIT NONE
  PUBLIC
  SAVE

  INTEGER, PARAMETER :: N = 70 !number of forward/adjoint model equations
  INTEGER, PARAMETER :: P = 2 !number of parameters
  DOUBLE PRECISION :: p1 !1st parameter (diffusion coefficient)
  DOUBLE PRECISION :: p2 !2nd parameter (convection coefficient)

  !initial and final integration time moments
  DOUBLE PRECISION :: t0, tF
  !spatial grid boundaries
  DOUBLE PRECISION, PARAMETER :: xStart = 0.d0, xEnd = 2.d0
  !mesh size
  DOUBLE PRECISION, PARAMETER :: dx = (xEnd - xStart) / (N+1)

  !~~~> reference solution of the forward convection-diffusion PDE
  DOUBLE PRECISION :: Yref(N)

END MODULE CDIFF_PARAMS


!Optimization routine parameters
MODULE LBFGS_PARAMS

  IMPLICIT NONE
  PUBLIC
  SAVE

!~~~> L-BFGS routine parameters
! For more information on the parameters, see the help comments in setulb.f

! nmax (integer) : the dimension of the largest problem to be solved.
! mmax (integer) : the maximum number of limited memory corrections.
  INTEGER, PARAMETER :: nmax = 2
  INTEGER, PARAMETER :: mmax = 10

  CHARACTER*60       ::    task, csave
  LOGICAL            ::    lsave(4)
  INTEGER            ::     m, iprint, nbd(nmax), iwa(3*nmax), isave(44)
  DOUBLE PRECISION   ::     f, factr, pgtol, params(nmax), &
                       &    l(nmax), u(nmax), g(nmax), dsave(29), &
                       &    wa(2*mmax*nmax+4*nmax+12*mmax*mmax+12*mmax)

! maximum number of function/gradient evaluations
  INTEGER, PARAMETER :: MaxFGeval = 200

END MODULE LBFGS_PARAMS


PROGRAM MAIN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Driver for the convection-diffusion problem
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  USE PERTURBATION
  USE CDIFF_PARAMS
  USE LBFGS_PARAMS
  IMPLICIT NONE
  INTRINSIC EXP

  !integrator parameters
  DOUBLE PRECISION :: rstatus(20), rcntrl(20)
  INTEGER :: ierr, istatus(20), icntrl(20)
  DOUBLE PRECISION :: ad_rstatus(20), ad_rcntrl(20)
  INTEGER :: ad_ierr, ad_istatus(20), ad_icntrl(20)

  !variables
  !Y (N by 1 array, double precision) : solution of forward conv-diff PDE
  DOUBLE PRECISION :: Y(N)
  !ADY (N by 1 array, double precision) : solution of adjoint conv-diff PDE
  DOUBLE PRECISION :: ADY(N)
  !W (P by 1 array, double precision) : quadrature system solution
  DOUBLE PRECISION :: W(P)
  !G (double precision) : cost function value
  DOUBLE PRECISION :: Gval

  !helper array to setup initial conditions
  DOUBLE PRECISION :: x(N)
  INTEGER :: i


  INTERFACE

      SUBROUTINE GRADEVAL(N, P, Y, ADY, W, TIN, TOUT, &
            ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U, &
            adICNTRL_U, adRCNTRL_U, adISTATUS_U, adRSTATUS_U, adIERR_U)

          INTEGER :: N
          INTEGER :: P
          DOUBLE PRECISION, INTENT(INOUT) :: Y(N), ADY(N), W(P)
          DOUBLE PRECISION, INTENT(INOUT) :: TIN  ! Start Time
          DOUBLE PRECISION, INTENT(INOUT) :: TOUT ! End Time
          ! Optional input parameters and statistics
          INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
          DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
          INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
          DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
          INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U

          INTEGER,       INTENT(IN),  OPTIONAL :: adICNTRL_U(20)
          DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: adRCNTRL_U(20)
          INTEGER,       INTENT(OUT), OPTIONAL :: adISTATUS_U(20)
          DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: adRSTATUS_U(20)
          INTEGER,       INTENT(OUT), OPTIONAL :: adIERR_U

      END SUBROUTINE GRADEVAL

      SUBROUTINE GET_YREF(N, Yref, TIN, TOUT, &
                  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U)

          INTEGER :: N
          DOUBLE PRECISION, INTENT(INOUT) :: Yref(N)
          DOUBLE PRECISION, INTENT(INOUT) :: TIN  ! Start Time
          DOUBLE PRECISION, INTENT(INOUT) :: TOUT ! End Time
          ! Optional input parameters and statistics
          INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
          DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
          INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
          DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
          INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U

      END SUBROUTINE GET_YREF

  END INTERFACE

  !set integrator tuning parameters to zero
  icntrl(:) = 0
  ad_icntrl(:) = 0
  rcntrl(:) = 0.d0
  ad_rcntrl(:) = 0.d0

  !set initial and final integration time
  t0 = 0.d0
  tF = 1.d0

  !set initial conditions for Y and Yref
  DO i=1,N
      x(i) = xStart + dx*i
      Y(i) = x(i)*(2.d0 - x(i))*exp(2.d0*x(i))
      Yref(i) = x(i)*(2.d0 - x(i))*exp(2.d0*x(i))
  END DO

  !initial condition for the adjoint PDE
  W(1:P) = 0.d0

!~~~> Set the nominal values for the parameters
  p1 = 1.d0
  p2 = 0.5d0
!~~~> Compute the reference solution Yref
  CALL GET_YREF(N, Yref, t0, tF, icntrl, rcntrl, istatus, rstatus, ierr)

  IF (ierr .ne. 1) THEN 
       PRINT *, 'Error occured during reference solution computation. Exiting.'
       STOP 0
  END IF

! !~~~> Code to estimate the gradient of G(tF) via finite differences
! !~~~> Uncomment to use
!   p1 = 1.5d0
!   p2 = 1.d0
!   !~~~> compute the gradient w
!   CALL GRADEVAL(N, P, Y, ADY, W, t0, tF, &
!               icntrl, rcntrl, istatus, rstatus, ierr, &
!               ad_icntrl, ad_rcntrl, ad_istatus, ad_rstatus, ad_ierr)
!   PRINT *, 'The gradient @tF: ', W
!   CALL GEVAL(Y,Gval)
!   PRINT *, 'Cost function : ', Gval
! 
!   p1 = 1.5d0 + dp1
!   p2 = 1.d0
!   DO i=1,N
!       Y(i) = x(i)*(2.d0 - x(i))*exp(2.d0*x(i))
!   END DO
!   CALL GET_YREF(N, Y, t0, tF, icntrl, rcntrl, istatus, rstatus, ierr)
!   CALL GEVAL(Y,Gval)
!   PRINT *, 'Cost function : ', Gval
! 
!   p1 = 1.5d0
!   p2 = 1.d0 + dp2
!   DO i=1,N
!       Y(i) = x(i)*(2.d0 - x(i))*exp(2.d0*x(i))
!   END DO
!   CALL GET_YREF(N, Y, t0, tF, icntrl, rcntrl, istatus, rstatus, ierr)
!   CALL GEVAL(Y,Gval)
!   PRINT *, 'Cost function : ', Gval


!~~~> Optimization problem:

!~~~> Set initial values for the parameters
  p1 = 3d0
  p2 = 3d0

!~~~> output f and ||proj g|| at every iteration
  iprint = 1
!~~~> tolerances in the stopping criteria
  factr = 1e+4
  pgtol = 1e-10

!~~~> p = 2 is the dimension of our optimization problem (2 parameters)
!~~~> m = number of limited memory corrections stored (must have m <= mmax)
  m=6

!~~~> Both variables are constrained
!~~~> We set both lower and upper bounds
  nbd(1)=2
  nbd(2)=2

!~~~> We impose 0.01 <= p1 <= 5
!           0.01 <= p2 <= 5
  l(1) = 0.01d0
  l(2) = 0.01d0
  u(1) = 5.d0
  u(2) = 5.d0

!~~~> We now define the starting point
  params(1) = p1
  params(2) = p2

!~~~> We start the iteration by initializing task = 'START'.
  task = 'START'

 222  CONTINUE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The call to the L-BFGS-B optimization subroutine SETULB
! We use the L-BFGS-B implementation by Ciyou Zhu,
! Richard Byrd and Jorge Nocedal.:
!
!   subroutine setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa,
!     +                 task, iprint, csave, lsave, isave, dsave)
!
!
!     NEOS, November 1994. (Latest revision April 1997.)
!     Optimization Technology Center.
!     Argonne National Laboratory and Northwestern University.
!     Written by
!                        Ciyou Zhu
!     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
!     References:
!
!       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
!       memory algorithm for bound constrained optimization'',
!       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!
!       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
!       limited memory FORTRAN code for solving bound constrained
!       optimization problems'', Tech. Report, NAM-11, EECS Department,
!       Northwestern University, 1994.
!
!                           *  *  *
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  CALL setulb(p,m,params,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint, &
            & csave,lsave,isave,dsave)

  !~~~> setulb requests the value of f and the gradient g at the current point
  ! (params(1), params(2))
  IF (task(1:2) .eq. 'FG') THEN

      !~~~> reinitialize Y
      DO i=1,N
         Y(i) = x(i)*(2.d0 - x(i))*exp(2.d0*x(i))
      END DO
      !~~~> reinitialize the gradient W
      W = 0.d0
      !~~~> reset the parameters p1 and p2
      p1 = params(1)
      p2 = params(2)

      !~~~> compute the gradient w
      CALL GRADEVAL(N, P, Y, ADY, W, t0, tF,   &
                    icntrl, rcntrl, istatus,   &
                    rstatus, ierr,             &
                    ad_icntrl, ad_rcntrl,      &
                    ad_istatus, ad_rstatus, ad_ierr)

      IF ((ierr .ne. 1) .or. (ad_ierr .ne. 1)) THEN 
            PRINT *, 'Error occured during gradient computation. Exiting.'
            STOP 0
      END IF

      !~~~> copy the gradient w into g
      g(1:2) = W(1:2)
      !~~~> evaluate the cost function and copy value in f
      CALL GEVAL(Y,f)

      !~~~> go back to the minimization subroutine
      GOTO 222

  ENDIF

  !setulb has returned with a new iterate
  IF (task(1:5) .eq. 'NEW_X')  THEN

     !~~~> reset the parameters p1 and p2
     p1 = params(1)
     p2 = params(2)

     !~~~> terminate the run if the total number of function/gradient
     ! evaluation exceeds the threshold value MaxFGeval
     IF (isave(34) .GE. MaxFGeval) THEN

         task='STOP: REACHED MAXIMUM NUMBER OF f/g EVALUATIONS'

     END IF

     IF (dsave(13) .LE. 1.d-8*(1.0d0 + abs(f))) THEN

         task = 'STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL.'

     END IF

     IF (task(1:4) .eq. 'STOP') then
         PRINT *, task
         PRINT *, 'Final parameters: '
         WRITE (6,'((1x,1p, 6(1x,d10.4)))') (params(i),i = 1,2)
     END IF

     !~~~> go back to the minimization subroutine
     GOTO 222

  END IF

!~~~>  If task is neither FG nor NEW_X we terminate execution
!     and exit the program.

  PRINT *, 'Final parameters:'
  WRITE (6,'((1x,1p, 6(1x,d10.4)))') (params(i),i = 1,2)

END PROGRAM MAIN



SUBROUTINE GEVAL(Y,Gval)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Evaluates the cost function
!      G(t^F) = 0.5 * \int_{xStart}^{xEnd} (y - y_ref)**2 dx
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  USE CDIFF_PARAMS
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: Y(N)
  DOUBLE PRECISION, INTENT(OUT) :: Gval
  INTEGER :: i

  Gval = 0.5d0 * (Y(1)-Yref(1)) * (Y(1)-Yref(1)) &
        + 0.5d0 * (Y(N)-Yref(N)) * (Y(N)-Yref(N))

  DO i=2,N-1
     Gval = Gval + (Y(i)-Yref(i)) * (Y(i)-Yref(i))
  END DO

  Gval = 0.5 * Gval * dx

END SUBROUTINE 



SUBROUTINE GRADEVAL(N, P, Y, ADY, W, TIN, TOUT, &
      ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U, &
      adICNTRL_U, adRCNTRL_U, adISTATUS_U, adRSTATUS_U, adIERR_U)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Sample wrapper for RKINT and RKINT_ADJDR. Initializes the user-defined
!  integrator parameters, calls RKINT and RKINT_ADJDR and returns the
!  integration results (Y, ADY, W).
!
! Inputs:
!      N (integer) : size of the ODE system/forward model state
!      Y (N by 1 array, double precision) : ODE system state vector at TIN
!      TIN (double precision) : Initial integration time
!      TOUT (double precision) : Final integration time
!      ICNTRL_U (20 by 1 array, integer)
!      RCNTRL_U (20 by 1 array, double precision)
!      adICNTRL_U (20 by 1 array, integer)
!      adRCNTRL_U (20 by 1 array, double precision)
! Outputs
!      Y (N by 1 array, double precision) : ODE system state vector at TOUT
!      ISTATUS_U (20 by 1 array, integer)
!      RSTATUS_U (20 by 1 array, double precision)
!      adISTATUS_U (20 by 1 array, integer)
!      adRSTATUS_U (20 by 1 array, double precision)
!      IERR_U (integer) : forward model integrator error code on exit
!                    1 : integration completed successfully
!                    < 0 : an error occured. See the accompanying error
!                        message for further details.
!      adIERR_U (integer) : adjoint model integrator error code on exit
!                    1 : integration completed successfully
!                    < 0 : an error occured. See the accompanying error
!                        message for further details.
!
!   For more information on ICNTRL_U, RCNTRL_U, ISTATUS_U and RSTATUS_U,
!   please see the source code for RKINT in ../src/erk.f90
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Use the modules containing the forward and adjoint Runge-Kutta solvers
   USE RK_MODULE
   USE RKADJDR_MODULE
! Use the memory buffers
   USE TAPES
! Use the convection-diffusion problem parameters
   USE CDIFF_PARAMS, ONLY : dx, Yref
   IMPLICIT NONE
! External subroutines called by forward/adjoint model integrators
   EXTERNAL FEVAL, ADJRHS, QUADRHS, JACEVAL_Y, DY0DP

   INTEGER :: N
   INTEGER :: P
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N), ADY(N), W(P)
   DOUBLE PRECISION, INTENT(INOUT) :: TIN  ! Start Time
   DOUBLE PRECISION, INTENT(INOUT) :: TOUT ! End Time
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U

   INTEGER,       INTENT(IN),  OPTIONAL :: adICNTRL_U(20)
   DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: adRCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: adISTATUS_U(20)
   DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: adRSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: adIERR_U


   DOUBLE PRECISION :: RCNTRL(20), RSTATUS(20)
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR

   DOUBLE PRECISION :: adRCNTRL(20), adRSTATUS(20)
   INTEGER       :: adICNTRL(20), adISTATUS(20), adIERR

   DOUBLE PRECISION :: ATOL(N), RTOL(N)
   DOUBLE PRECISION :: STEPMIN
   DOUBLE PRECISION :: adATOL(N+P), adRTOL(N+P)
   DOUBLE PRECISION :: adSTEPMIN

   DOUBLE PRECISION :: T1, T2 !for timing purposes only
   LOGICAL :: CheckptEnabled

   INTEGER :: Nd, Nc, TapeSize, MaxFwdStps, MaxAdjStps

   INTEGER, PARAMETER :: RK5 = 4, &  !DOPRI5(4)
                         RK2 = 1, &  !Fehlberg RK2(3)
                         RK3 = 2, &  !RK3(2)
                         RK4 = 3, &  !RK4(3)
                         RK6 = 5, &  !RK6(5)
                         RK8 = 6     !DOPRI8(6)

   !Checkpointing will be used
   CheckptEnabled = .TRUE.
   Nd = 1000 !checkpoint every Nd steps

   !Maximum number of forward/adjoint integration steps
   MaxFwdStps = 50000
   MaxAdjStps = 50000

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.d0
   ISTATUS(:) = 0
   RSTATUS(:) = 0.d0
   adICNTRL(:)  = 0
   adRCNTRL(:)  = 0.d0
   adISTATUS(:) = 0
   adRSTATUS(:) = 0.d0

   !~~~> fine-tune the integrator:
   ICNTRL(1) = 0    ! 0 - vector tolerances, 1 - scalars
   adICNTRL(1) = 0

   ICNTRL(2) = MaxFwdStps
   adICNTRL(2) = MaxAdjStps

   !Choose a Runge-Kutta method
   ICNTRL(3) = RK5
   adICNTRL(3) = RK5

!~~~> Default value: adICNTRL(4) = 0: dense output enabled
!~~~> Uncomment any of the following two lines to enable Hermite interpolation
   !adICNTRL(4) = 1 !3rd order Hermite
   !adICNTRL(4) = 2 !5th order Hermite

!~~~> Problem is parameter-dependent (eps)
   adICNTRL(5) = 1

!~~~> Enable checkpointing and memory buffering
   ICNTRL(4) = 1 !enable file checkpointing
   ICNTRL(5) = 1 !enable memory buffering
   adICNTRL(6) = 1 !inform RKINT_ADJDR that checkpointing
                   !was used in FWD mode

   ! --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
   RTOL(:)=1.0D-12
   ATOL(:)=RTOL(:)
   adRTOL(:)=1.0D-10
   adATOL(:)=adRTOL(:)

   ! If optional parameters are given, and if they are >0,
   ! then they overwrite default settings.
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF

   IF (PRESENT(adICNTRL_U)) THEN
     WHERE(adICNTRL_U(:) > 0) adICNTRL(:) = adICNTRL_U(:)
   END IF
   IF (PRESENT(adRCNTRL_U)) THEN
     WHERE(adRCNTRL_U(:) > 0) adRCNTRL(:) = adRCNTRL_U(:)
   END IF

   IF (CheckptEnabled) THEN
       TapeSize = 2*Nd
       CALL tapes_Allocate(CheckptEnabled,TapeSize,N,16)
   ELSE
       TapeSize = ICNTRL(2)
       CALL tapes_Allocate(CheckptEnabled,TapeSize,N,16)
   END IF

   !Start clock
   CALL cpu_time(T1)
   !Forward integration
   CALL RKINT(N,Y,FEVAL,TIN,TOUT,  &
         Nd, Nc,                   &
         ATOL,RTOL,                &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
   !Stop clock
   CALL cpu_time(T2)

!    PRINT *,'Forward model integration took ', (T2-T1)*1000.d0 , ' ms.'

   !~~~> If optional parameters are given for output they
   ! are updated with the return information
   !~~~> Update statistics for FWD integration and reset ISTATUS
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   ISTATUS(:) = 0

   !~~~> Check for errors during the adjoint model integration
   IF (IERR .NE. 1) THEN
       PRINT *, 'RKINT returned with an error!'
       CALL rk_DeallocateTapes
       RETURN
   END IF

!    IF (CheckptEnabled) THEN
!        PRINT *, 'Wrote ', Nc, ' checkpoint(s) during FWD integration'
!    END IF

   !Initialize the adjoint variables (at TOUT)
   !We have g(y) = 0.5 * (y - y_ref)**2, therefore ADY(tF) = dg/dy = y - y_ref
   ADY = Y - Yref

   !Start clock
   CALL cpu_time(T1)
   !Adjoint integration
   CALL RKINT_ADJDR(N, ADY, ADJRHS, JACEVAL_Y,               &
           P, W, QUADRHS, DY0DP, FEVAL,                      &
           TIN, TOUT, Nc,                                    &
           ATOL, RTOL,                                       &
           adATOL,adRTOL,                                    &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR,               &
           adRCNTRL,adICNTRL,adRSTATUS,adISTATUS,adIERR)
   !Stop clock
   CALL cpu_time(T2)

!  PRINT *,'Adjoint model integration took ', (T2-T1)*1000.d0 , ' ms.'

   !~~~> Check for errors during the adjoint model integration
   IF (adIERR .NE. 1) THEN
       PRINT *, 'RKINT_ADJ returned with an error!'
       CALL rk_DeallocateTapes
       RETURN
   END IF

   !~~~> Deallocate the memory buffers
   CALL rk_DeallocateTapes

   ! If optional parameters are given for output they
   ! are updated with the return information
   ! add statistics for FWD integrations on subintervals (if needed)
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS_U(:) + ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

   IF (PRESENT(adISTATUS_U)) adISTATUS_U(:) = adISTATUS(:)
   IF (PRESENT(adRSTATUS_U)) adRSTATUS_U(:) = adRSTATUS(:)
   IF (PRESENT(adIERR_U))    adIERR_U       = adIERR

   STEPMIN = RSTATUS(Nhexit)
   adSTEPMIN = RSTATUS(adNhexit)

END SUBROUTINE GRADEVAL



SUBROUTINE GET_YREF(N, Yref, TIN, TOUT,            &
                    ICNTRL_U, RCNTRL_U, ISTATUS_U, &
                    RSTATUS_U, IERR_U)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Computes the reference solution Yref for p1 = 1.d0 and p2 = 0.5d0
!
! Inputs:
!      N (integer) : size of the ODE system/forward model state
!      Yref (N by 1 array, double precision) : ODE system state vector at TIN
!      TIN (double precision) : Initial integration time
!      TOUT (double precision) : Final integration time
!      ICNTRL_U (20 by 1 array, integer)
!      RCNTRL_U (20 by 1 array, double precision)
! Outputs
!      Yref (N by 1 array, double precision) : ODE system state vector at TOUT
!      ISTATUS_U (20 by 1 array, integer)
!      RSTATUS_U (20 by 1 array, double precision)
!      IERR_U (integer) : error code on exit
!                    1 : integration completed successfully
!                    < 0 : an error occured. See the accompanying error 
!                        message for further details.
!
!   For more information on ICNTRL_U, RCNTRL_U, ISTATUS_U and RSTATUS_U,
!   please see the source code for RKINT in ../src/erk.f90
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   USE RK_MODULE
   USE RKADJDR_MODULE
   USE TAPES
   IMPLICIT NONE
   EXTERNAL FEVAL

   INTEGER :: N
   DOUBLE PRECISION, INTENT(INOUT) :: Yref(N)
   DOUBLE PRECISION, INTENT(INOUT) :: TIN  ! Start Time
   DOUBLE PRECISION, INTENT(INOUT) :: TOUT ! End Time
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U

   DOUBLE PRECISION :: RCNTRL(20), RSTATUS(20)
   INTEGER       :: ICNTRL(20), ISTATUS(20), IERR

   DOUBLE PRECISION :: ATOL(N), RTOL(N)
   DOUBLE PRECISION :: STEPMIN

   INTEGER, PARAMETER :: Nd = 0
   INTEGER :: Nc
   INTEGER, PARAMETER :: MaxFwdSteps = 50000
   DOUBLE PRECISION :: t1, t2 !for timing purposes

   INTEGER, PARAMETER :: RK5 = 4, &  !DOPRI5(4)
                         RK2 = 1, &  !Fehlberg RK2(3)
                         RK3 = 2, &  !RK3(2)
                         RK4 = 3, &  !RK4(3)
                         RK6 = 5, &  !RK6(5)
                         RK8 = 6     !DOPRI8(6)

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.d0
   ISTATUS(:) = 0
   RSTATUS(:) = 0.d0

   !~~~> fine-tune the integrator:
   ICNTRL(1) = 0    ! 0 - vector tolerances, 1 - scalars
   ICNTRL(2) = MaxFwdSteps
   ICNTRL(3) = RK5

! --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
! --- Compute Yref to high accuracy
   RTOL(:)=1.0D-12
   ATOL(:)=RTOL(:)

!~~~> If optional parameters are given, and if they are >0,
! then they overwrite default settings.
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF

   !Start clock
   CALL cpu_time(t1)

   !Forward PDE integration
   CALL RKINT(N,Yref,FEVAL,TIN,TOUT,           &
         Nd, Nc,                               &
         ATOL,RTOL,                            &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)

   !Stop clock
   CALL cpu_time(t2)
   PRINT *,'Computing the reference value Yref took ', (T2-T1)*1000.d0 , ' ms.'

   STEPMIN = RSTATUS(Nhexit)
   !~~~> If optional parameters are given for output they 
   ! are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U)) IERR_U = IERR

END SUBROUTINE GET_YREF



SUBROUTINE FEVAL(N, T1, Y1, F1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! FEVAL computes the right-hand side of the forward ODE Y' = F(T,Y)
! at (T1,Y1), i.e. F1 = F(T1,Y1).
!
! Inputs:
!      N (integer) : size of the ODE system
!      T1 (double precision) : current integration time
!      Y1 (N by 1 array, double precision) : ODE numerical solution at time T1
!
! Output:
!      F1 (N by 1 array, double precision) : ODE right-hand side at time T1
!                                            F1 = F(T1,Y1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  USE CDIFF_PARAMS, ONLY : dx, p1, p2
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, INTENT(IN) :: T1, Y1(N)
  DOUBLE PRECISION, INTENT(OUT) :: F1(N)

  INTEGER :: i

  F1(1) = p1*(Y1(2) - 2.d0*Y1(1))/(dx**2) + &
         & p2*Y1(2)/(2.d0*dx)

  DO i = 2,N-1
     F1(i) =  p1*(Y1(i+1) - 2.d0*Y1(i) + Y1(i-1))/(dx**2) + &
            & p2*(Y1(i+1) - Y1(i-1))/(2.d0*dx)
  END DO

  F1(N) = p1*(-2.d0*Y1(N) + Y1(N-1))/(dx**2) + &
            & p2*(-Y1(N-1))/(2.d0*dx)

END SUBROUTINE FEVAL



SUBROUTINE ADJRHS(N, T1, Y1, ADY1, ADJ1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Evaluates the RHS of the adjoint ODE (the spatial discretization of the
! adjoint PDE) at (T1,Y1,ADY1)
!
! Inputs:
!      N (integer) : size of the ODE system
!      T1 (double precision) : current integration time
!      Y1 (N by 1 array, double precision) : forward ODE numerical solution at time T1
!    ADY1 (N by 1 array, double precision) : adjoint ODE numerical solution at time T1
! Output:
!    ADJ1 (N by 1 array, double precision) : adjoint ODE right-hand side at time T1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  USE CDIFF_PARAMS, ONLY : dx, p1, p2
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, INTENT(IN) :: T1, Y1(N), ADY1(N)
  DOUBLE PRECISION, INTENT(OUT) :: ADJ1(N)

  INTEGER :: i

  ADJ1(1) = -p1*(ADY1(2) - 2.d0*ADY1(1))/(dx**2) + &
         & p2*ADY1(2)/(2.d0*dx)

  DO i = 2,N-1
      ADJ1(i) = -p1*(ADY1(i+1) - 2.d0*ADY1(i) + ADY1(i-1))/(dx**2) + &
            & p2*(ADY1(i+1) - ADY1(i-1))/(2.d0*dx)
  END DO

  ADJ1(N) = -p1*(-2.d0*ADY1(N) + ADY1(N-1))/(dx**2) + &
            & p2*(-ADY1(N-1))/(2.d0*dx)

END SUBROUTINE ADJRHS



SUBROUTINE JACEVAL_Y(N, T1, Y1, g_Y1, g_F1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Evaluates the RHS of the quadrature ODE system.
!
! Inputs:
!      N (integer) : size of the forward/adjoint ODE system
!      P (integer) : size of the quadrature ODE system
!      T1 (double precision) : current integration time
!      Y1 (N by 1 array, double precision) : forward ODE numerical solution at time T1
!    g_Y1 (N by 1 array, double precision) : TLM ODE numerical solution at time T1
! Output:
!    g_F1 (P by 1 array, double precision) : Jacobian-vector product (dF/dy)*g_Y1 at (T1,Y1)
!                                          This is useful when performing Hermite interpolation
!                                          of the forward model trajectory.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE CDIFF_PARAMS, ONLY : dx, p1, p2
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, INTENT(IN) :: T1, Y1(N), g_Y1(N)
  DOUBLE PRECISION, INTENT(OUT) :: g_F1(N)

  INTEGER :: i

  g_f1(1) = g_y1(2)*(p1/dx**2+p2/(2.d0*dx))+g_y1(1)*((-2)*p1/dx**2)
  DO i = 2, n-1
    g_f1(i) = g_y1(i-1)*(p1/dx**2-p2/(2.d0*dx))+g_y1(i+1)*(p1/dx**2+p2/(2.d0*dx))+g_y1(i)*((-2)*p1/dx**2)
  END DO
  g_f1(n) = g_y1(n-1)*(p1/dx**2-p2/(2.d0*dx))+g_y1(n)*((-2)*p1/dx**2)


END SUBROUTINE JACEVAL_Y



SUBROUTINE QUADRHS(N,P,T1,Y1,ADY1,W1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Evaluates the RHS of the quadrature ODE system.
!
! Inputs:
!      N (integer) : size of the forward/adjoint ODE system
!      P (integer) : size of the quadrature ODE system
!      T1 (double precision) : current integration time
!      Y1 (N by 1 array, double precision) : forward ODE numerical solution at time T1
!    ADY1 (N by 1 array, double precision) : adjoint ODE numerical solution at time T1
! Output:
!    W1 (P by 1 array, double precision) : quadrature ODE right-hand side at time T1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE CDIFF_PARAMS, ONLY : dx, Yref
  USE UTILS
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N, P
  DOUBLE PRECISION, INTENT(IN) :: T1
  DOUBLE PRECISION, INTENT(IN) :: Y1(N), ADY1(N)
  DOUBLE PRECISION, INTENT(OUT) :: W1(P)

  INTEGER :: i
  DOUBLE PRECISION :: tmp(N)

  tmp(1) = (Y1(2) - 2.d0*Y1(1))/dx
  DO i=2,N-1
     tmp(i) = (Y1(i+1) - 2.d0*Y1(i) + Y1(i-1))/dx
  END DO
  tmp(N) = (-2.d0*Y1(N) + Y1(N-1))/dx

!~~~> compute the first spatial integral using the trapezoidal rule
  W1(1) = 0.5d0*(tmp(1)*ADY1(1)*(Y1(1)-Yref(1)) + tmp(N)*ADY1(N))
  DO i=2,N-1
     W1(1) = W1(1) + tmp(i)*ADY1(i)
  END DO

  tmp(1) = 0.5d0*Y1(2)
  DO i=2,N-1
     tmp(i) = 0.5d0*(Y1(i+1) - Y1(i-1))
  END DO
  tmp(N) = -0.5d0*Y1(N-1)

!~~~> compute the second spatial integral using the trapezoidal rule
  W1(2) = 0.5d0*(tmp(1)*ADY1(1) + tmp(N)*ADY1(N))
  DO i=2,N-1
      W1(2) = W1(2) + tmp(i)*ADY1(i)
  END DO

  !CALL WSCAL(P,-1.d0,W1,1)
  W1 = -W1

END SUBROUTINE QUADRHS



SUBROUTINE DY0DP(N,P,T0,ADY0,GRAD0)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 !Computes the final update at T0 to the gradient dG/dp
 !i.e. (lambda(t0)*s(t0))^T where s(t0)=dy0/dp (the sensitivity of the
 !initial solution w.r.t. the P parameters - a N x P matrix)
 !
 !The initial conditions are parameter-independent, so no correction
 !is needed here (the subroutine can simply return a zero vector).
 !
 !Inputs:
 !  N = size of the forward ODE system
 !  P = number of parameters in the system
 !  ADY0 = the adjoint solution at the initial time
 !Output:
 !  GRAD0 = Correction to the gradient dG/dp(t0); needed if the initial
 !conditions of the forward ODE system are parameter-dependent
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N, P
  DOUBLE PRECISION, INTENT(IN) :: T0
  DOUBLE PRECISION, INTENT(IN)  :: ADY0(N)
  DOUBLE PRECISION, INTENT(OUT) :: GRAD0(P)

  GRAD0 = 0.d0

END SUBROUTINE DY0DP
