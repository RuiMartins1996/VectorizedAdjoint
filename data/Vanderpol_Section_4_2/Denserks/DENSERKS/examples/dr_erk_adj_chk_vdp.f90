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


PROGRAM MAIN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Driver for the first order adjoint of the scaled Van der Pol system
! eps = 0.001d0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  INTEGER, PARAMETER :: N = 2
  INTEGER, PARAMETER :: P = 1 !consider eps to be a parameter
  DOUBLE PRECISION :: t0, tF
  DOUBLE PRECISION :: rstatus(20), rcntrl(20)
  INTEGER :: ierr, istatus(20), icntrl(20)
  DOUBLE PRECISION :: ad_rstatus(20), ad_rcntrl(20)
  INTEGER :: ad_ierr, ad_istatus(20), ad_icntrl(20)
  DOUBLE PRECISION :: Y(N), ADY(N), W(P)
  DOUBLE PRECISION, PARAMETER :: eps = 0.001d0


  DOUBLE PRECISION, PARAMETER :: h = 1e-5

  INTERFACE

      SUBROUTINE INTEGRATE(N, P, Y, ADY, W, TIN, TOUT, &
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

      END SUBROUTINE INTEGRATE

  END INTERFACE

! Set the initial and final integration time
  t0 = 0.d0
  tF = 5e-1  

! Set initial conditions
  Y(1) = 2.d0
  Y(2) = -2.0d0 / 3.0d0 + 10.0d0 / (81.0 / eps) - 292.0d0 / (2187.0 / eps**2)

  ADY(1) = 1.d0
  ADY(2) = 0.d0

  W = 0.d0

! Set all integrator parameters to zero
  rcntrl(:) = 0.d0
  ad_rcntrl(:) = 0.d0
  icntrl(:) = 0
  ad_icntrl(:) = 0

  CALL INTEGRATE(N, P, Y, ADY, W, t0, tF, &
        & icntrl, rcntrl, istatus, rstatus, ierr, &
        & ad_icntrl, ad_rcntrl, ad_istatus, ad_rstatus, ad_ierr)

! --- PRINT FINAL SOLUTION
  PRINT *,'Final solution:'
  PRINT *,'Y(1) = ',  Y(1)
  PRINT *,'Y(2) = ',  Y(2)

  PRINT *,'No. of function calls: ', istatus(1)
  PRINT *,'No. of steps (FWD): '   , istatus(2)
  PRINT *,'No. of accepted steps: ', istatus(3)
  PRINT *,'No. of rejected steps: ', istatus(4)

  PRINT *,'*****************'
! --- PRINT ADJOINT SOLUTION
  PRINT *,'Adjoint solution:'
  PRINT *,'ADY(1) = ',  ADY(1)
  PRINT *,'ADY(2) = ',  ADY(2)
!--- PRINT STATISTICS
  PRINT *,'No. of Jacobian-(transpose)-vector products: ', ad_istatus(1)
  PRINT *,'No. of steps (ADJ): ',              ad_istatus(2)
  PRINT *,'No. of accepted steps: ',           ad_istatus(3)
  PRINT *,'No. of rejected steps: ',           ad_istatus(4)
  PRINT *,'No. of second derivative evaluations (Hermite):', ad_istatus(5)
!
  PRINT *,'Gradient d_yF/d_eps:'
  PRINT *,'W(1) = ',  W(1)

END PROGRAM MAIN


SUBROUTINE INTEGRATE(N, P, Y, ADY, W, TIN, TOUT, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U, &
  adICNTRL_U, adRCNTRL_U, adISTATUS_U, adRSTATUS_U, adIERR_U)

   USE RK_MODULE
   USE RKADJDR_MODULE
   IMPLICIT NONE
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
   ICNTRL(3) = RK8
   adICNTRL(3) = RK8

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
   adRTOL(:)=1.0D-12
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
       TapeSize = Nd + 1
       CALL rk_AllocateTapes(CheckptEnabled,TapeSize,N)
   ELSE
       TapeSize = ICNTRL(2)
       CALL rk_AllocateTapes(CheckptEnabled,TapeSize,N)
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
   PRINT *,'Forward model integration took ', (T2-T1)*1000.d0 , ' ms.'

   ! if optional parameters are given for output they 
   ! are updated with the return information
   ! update statistics for FWD integration and reset ISTATUS
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   ISTATUS(:) = 0

   !check for errors during the adjoint model integration
   IF (IERR .NE. 1) THEN
       PRINT *, 'RKINT returned with an error!'
       CALL rk_DeallocateTapes
       RETURN
   END IF

   IF (CheckptEnabled) THEN
       PRINT *, 'Wrote ', Nc, ' checkpoint(s) during FWD integration'
   END IF

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
   PRINT *,'Adjoint model integration took ', (T2-T1)*1000.d0 , ' ms.'

   !check for errors during the adjoint model integration
   IF (adIERR .NE. 1) THEN
       PRINT *, 'RKINT_ADJ returned with an error!'
       CALL rk_DeallocateTapes
       RETURN
   END IF

   CALL rk_DeallocateTapes

   STEPMIN = RSTATUS(Nhexit)
   ! if optional parameters are given for output they 
   ! are updated with the return information
   ! add statistics for FWD integrations on subintervals (if needed)
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS_U(:) + ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

   adSTEPMIN = RSTATUS(adNhexit)
   ! if optional parameters are given for output they 
   ! are updated with the return information
   IF (PRESENT(adISTATUS_U)) adISTATUS_U(:) = adISTATUS(:)
   IF (PRESENT(adRSTATUS_U)) adRSTATUS_U(:) = adRSTATUS(:)
   IF (PRESENT(adIERR_U))    adIERR_U       = adIERR

END SUBROUTINE INTEGRATE



SUBROUTINE FEVAL(N,T1,Y1,F1)
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

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, INTENT(IN) :: T1, Y1(N)
  DOUBLE PRECISION, INTENT(OUT) :: F1(N)
  DOUBLE PRECISION, PARAMETER :: eps = 0.001d0

  !we're solving Van der Pol's equation
  !as a variation, the user can consider eps to be a parameter and compute
  !dy(tF)/dp or dG/dp for a given cost functional G(y)

  F1(1) = Y1(2)
  F1(2) = ((1.d0 - Y1(1)**2)*Y1(2) - Y1(1)) / eps

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

  IMPLICIT NONE
  integer, intent(in) :: N
  double precision, intent(in) :: t1, y1(N), ady1(N)
  double precision, intent(out) :: adj1(N)
  double precision, parameter :: eps = 0.001d0

  adj1(1) = -(-2.d0*y1(1)*y1(2)-1.d0)*ady1(2)/eps
  adj1(2) = -ady1(1) - (1.d0-y1(1)**2)*ady1(2)/eps

END SUBROUTINE ADJRHS



SUBROUTINE JACEVAL_Y(N, T1, Y1, DY1, JAC1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Evaluates the RHS of the quadrature ODE system.
!
! Inputs:
!      N (integer) : size of the forward/adjoint ODE system
!      P (integer) : size of the quadrature ODE system
!      T1 (double precision) : current integration time
!      Y1 (N by 1 array, double precision) : forward ODE numerical solution at time T1
!     DY1 (N by 1 array, double precision) : TLM ODE numerical solution at time T1
! Output:
!    JAC1 (P by 1 array, double precision) : Jacobian-vector product (dF/dy)*DY1 at (T1,Y1)
!                                          This is useful when performing Hermite interpolation
!                                          of the forward model trajectory.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  IMPLICIT NONE
  integer, intent(in) :: N
  double precision, intent(in) :: t1, y1(N), dy1(N)
  double precision, intent(out) :: jac1(N)
  double precision, parameter :: eps = 0.001d0

  jac1(1) = dy1(2)
  jac1(2) = (-2.d0*y1(1)*y1(2)-1.d0)*dy1(1)/eps + (1.d0 - y1(1)**2)*dy1(2)/eps

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
  IMPLICIT NONE
  INTEGER :: N
  INTEGER :: P
  DOUBLE PRECISION, INTENT(IN) :: T1
  DOUBLE PRECISION, INTENT(IN) :: Y1(N), ADY1(N)
  DOUBLE PRECISION, INTENT(OUT) :: W1(P)
  DOUBLE PRECISION :: dfdp(2), eps

  !w' = -(df/dp)^T * w - (dg/dp)^T

  eps = 0.001d0
  dfdp(1) = 0.d0
  dfdp(2) = (-1.d0 / eps**2) * ((1.d0 - y1(1)**2)*y1(2) - y1(1))

  W1 = - dfdp(1)*ADY1(1) - dfdp(2)*ADY1(2)

END SUBROUTINE QUADRHS



SUBROUTINE DY0DP(N,P,T0,ADY0,GRAD0)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 !Computes the final update to the gradient dG/dp
 !i.e. (lambda(t0)*s(t0))^T where s(t0)=dy0/dp (the sensitivity of the
 !initial solution w.r.t. the P parameters - a N x P matrix)
 !If the initial conditions are parameter-independent, then no correction
 !is needed (the subroutine must simply return a zero vector).
 !
 ! Input:
 !  N = size of the forward ODE system
 !  P = number of parameters in the system
 !  T0 = initial integration time
 !  ADY0 = the adjoint solution at the initial time
 !
 ! Output:
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