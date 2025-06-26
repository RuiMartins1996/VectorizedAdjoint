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

MODULE PARAMETERS

  IMPLICIT NONE
  PUBLIC
  SAVE

  !np(integer) = number of control parameters
  INTEGER, PARAMETER :: np = 400
  !params (np by 1 array, double precision) = parameter values
  DOUBLE PRECISION :: params(np)
  !parameter perturbation vector
  DOUBLE PRECISION :: DP(np)

END MODULE PARAMETERS


PROGRAM MAIN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! First order sensitivity analysis for the parameter-dependent Van der Pol
! ODE system (for details, see section 5.2.2 in the DENSERKS paper).
! We compute two values for the gradient dg/dp at different points:
!
!              (dg/dp)(params)  and (dg/dp)(params + eps*DP)
!
! these values can be used to approximate the Hessian-vector product:
!
!     (d2g/dp2)*DP \approx (1/eps)*((dg/dp)(params + eps*DP) - (dg/dp)(params))
!
! The accuracy of this finite difference approximation depends on the choice of eps.
!
! The second order adjoint method allows us to obtain the value of the Hess-Vector
! product with controlled accuracy (via the user-selected integration tolerances).
!
! One can use a call to RKINT_SOADR_M to solve the two SOA systems more efficiently.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE PARAMETERS
  USE UTILS
  IMPLICIT NONE

  INTEGER, PARAMETER :: N = 3
  INTEGER, PARAMETER :: P = np
  DOUBLE PRECISION :: rstatus(20), rcntrl(20)
  INTEGER :: ierr, istatus(20), icntrl(20)
  DOUBLE PRECISION :: ad_rstatus(20), ad_rcntrl(20)
  INTEGER :: ad_ierr, ad_istatus(20), ad_icntrl(20)
  DOUBLE PRECISION :: Y(N), ADY(N)
  DOUBLE PRECISION :: W(P)
  !initial and final integration time
  DOUBLE PRECISION :: t0
  DOUBLE PRECISION :: tF

  !Finite differences "h"
  DOUBLE PRECISION, PARAMETER :: eps = 0.0001d0
  INTEGER :: i
  !FD approximation of the Hessian-Vector product
  DOUBLE PRECISION :: HessVFD(P)

  INTRINSIC SQRT, REAL

  INTERFACE 

    SUBROUTINE INTEGRATE(N, P, Y, ADY, W, TIN, TOUT,                 &
        ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U,            &
        adICNTRL_U, adRCNTRL_U, adISTATUS_U, adRSTATUS_U, adIERR_U)

           INTEGER :: N
           INTEGER :: P
           DOUBLE PRECISION, INTENT(INOUT) :: Y(N), ADY(N)
           DOUBLE PRECISION, INTENT(INOUT) :: W(P)
           DOUBLE PRECISION  :: TIN  ! Start Time
           DOUBLE PRECISION  :: TOUT ! End Time
           ! Optional input parameters and statistics
           INTEGER,        OPTIONAL :: ICNTRL_U(20)
           DOUBLE PRECISION,OPTIONAL :: RCNTRL_U(20)
           INTEGER,        OPTIONAL :: ISTATUS_U(20)
           DOUBLE PRECISION, OPTIONAL :: RSTATUS_U(20)
           INTEGER,        OPTIONAL :: IERR_U

           INTEGER,        OPTIONAL :: adICNTRL_U(20)
           DOUBLE PRECISION, OPTIONAL :: adRCNTRL_U(20)
           INTEGER,       OPTIONAL :: adISTATUS_U(20)
           DOUBLE PRECISION, OPTIONAL :: adRSTATUS_U(20)
           INTEGER,       OPTIONAL :: adIERR_U

    END SUBROUTINE INTEGRATE

  END INTERFACE 

  t0 = 0.d0
  tF = 5.d0

  !set initial conditions
  Y(1) = 0.d0
  Y(2) = 1.d0
  Y(3) = 0.d0

  !ADY(t0) = g_y(tF,y,p) = (0,0,1)^T
  ADY(1) = 0.d0
  ADY(2) = 0.d0
  ADY(3) = 1.d0

  CALL SET2ZERO(P,W)

  DO i = 1,P
      params(i) = 1.d0 / P
  END DO

  DO i=1,P
      DP(i) = 1.d0 / i
  END DO

! Set all integrator parameters to zero
  rcntrl(:) = 0.d0
  ad_rcntrl(:) = 0.d0
  icntrl(:) = 0
  ad_icntrl(:) = 0

  CALL INTEGRATE(N, P, Y, ADY, W, t0, tF,                         &
        & icntrl, rcntrl, istatus, rstatus, ierr,                 &
        & ad_icntrl, ad_rcntrl, ad_istatus, ad_rstatus, ad_ierr)

! Reinitialize Y, ADY, W and recompute the gradient and then validate
! the Hessian-vector product using finite differences

! --- PRINT FINAL SOLUTION
  PRINT *,'Final solution:'
  PRINT *,'Y(1) = ',  Y(1)
  PRINT *,'Y(2) = ',  Y(2)
  PRINT *,'Y(3) = ',  Y(3)

  PRINT *,'No. of function calls: ', istatus(1)
  PRINT *,'No. of steps (FWD): '   , istatus(2)
  PRINT *,'No. of accepted steps: ', istatus(3)
  PRINT *,'No. of rejected steps: ', istatus(4)

  PRINT *,'*****************'
! --- PRINT FIRST ORDER ADJOINT SOLUTION
  PRINT *,'First order adjoint solution:'
  PRINT *,'ADY(1) = ',  ADY(1)
  PRINT *,'ADY(2) = ',  ADY(2)
  PRINT *,'ADY(3) = ',  ADY(3)

! --- PRINT SECOND ORDER ADJOINT SOLUTION
!  PRINT *,'Second order adjoint solution:'
!  PRINT *,'AD2Y(1) = ',  AD2Y(1)
!  PRINT *,'AD2Y(2) = ',  AD2Y(2)
!  PRINT *,'AD2Y(3) = ',  AD2Y(3)
!  PRINT *, '****************'

!--- PRINT STATISTICS
  PRINT *,'No. of Jacobian-(transpose)-vector products: ', ad_istatus(1)
  PRINT *,'No. of steps (ADJ): ',              ad_istatus(2)
  PRINT *,'No. of accepted steps: ',           ad_istatus(3)
  PRINT *,'No. of rejected steps: ',           ad_istatus(4)
  PRINT *,'No. of second derivative evaluations (Hermite):', ad_istatus(5)

! PRINT *, 'W = '
! PRINT *, W

! Reset initial conditions for the forward model
  Y(1) = 0.d0
  Y(2) = 1.d0
  Y(3) = 0.d0

! Add the perturbation in the parameters to their nominal values
! params = params + eps*DP
  CALL WAXPY(P,eps,DP,1,params,1)

!Reinitialize the adjoint variable
  ADY(1) = 0.d0
  ADY(2) = 0.d0
  ADY(3) = 1.d0

  CALL SET2ZERO(P,HessVFD)

! Set all integrator parameters to zero
  rcntrl(:) = 0.d0
  ad_rcntrl(:) = 0.d0
  icntrl(:) = 0
  ad_icntrl(:) = 0

  CALL INTEGRATE(N, P, Y, ADY, HessVFD, t0, tF,                    &
          & icntrl, rcntrl, istatus, rstatus, ierr,                &
          & ad_icntrl, ad_rcntrl, ad_istatus, ad_rstatus, ad_ierr)

  PRINT *, 'Final solution'
  PRINT *,'Y(1) = ',  Y(1)
  PRINT *,'Y(2) = ',  Y(2)
  PRINT *,'Y(3) = ',  Y(3)

  PRINT *,'No. of function calls: ', istatus(1)
  PRINT *,'No. of steps (FWD): '   , istatus(2)
  PRINT *,'No. of accepted steps: ', istatus(3)
  PRINT *,'No. of rejected steps: ', istatus(4)

  PRINT *,'*****************'
! --- PRINT FIRST ORDER ADJOINT SOLUTION
  PRINT *,'First order adjoint solution:'
  PRINT *,'ADY(1) = ',  ADY(1)
  PRINT *,'ADY(2) = ',  ADY(2)
  PRINT *,'ADY(3) = ',  ADY(3)

! --- PRINT SECOND ORDER ADJOINT SOLUTION
!  PRINT *,'Second order adjoint solution:'
!  PRINT *,'AD2Y(1) = ',  AD2Y(1)
!  PRINT *,'AD2Y(2) = ',  AD2Y(2)
!  PRINT *,'AD2Y(3) = ',  AD2Y(3)
!  PRINT *, '****************'

!--- PRINT STATISTICS
  PRINT *,'No. of Jacobian-(transpose)-vector products: ', ad_istatus(1)
  PRINT *,'No. of steps (ADJ): ',              ad_istatus(2)
  PRINT *,'No. of accepted steps: ',           ad_istatus(3)
  PRINT *,'No. of rejected steps: ',           ad_istatus(4)
  PRINT *,'No. of second derivative evaluations (Hermite):', ad_istatus(5)

!  PRINT *, 'W = '
!  PRINT *, HessVFD

  CALL WAXPY(P,-1.d0,W,1,HessVFD,1)
  CALL WSCAL(P,1.d0/eps,HessVFD,1)

!   PRINT *, 'Finite difference approximation for the Hessian-vector ', &
!         'product HessVFD = (d2g/dp2)*DP :'
!   DO i=1,P
!       PRINT *, 'HessVFD(', i, ') = ', HessVFD(i)
!   END DO

END PROGRAM MAIN



SUBROUTINE INTEGRATE(N, P, Y, ADY, W, TIN, TOUT, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U, &
  adICNTRL_U, adRCNTRL_U, adISTATUS_U, adRSTATUS_U, adIERR_U)

!  Use the RKINT integrator for forward model integrations
   USE RK_MODULE
!  Use the RKINT_ADJDR wrapper to integrate the adjoint
   USE RKADJDR_MODULE
   IMPLICIT NONE
   EXTERNAL FEVAL, ADJRHS, JACV, DY0DP, QUADRHS

   INTEGER :: N
   INTEGER :: P
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N), ADY(N)
   DOUBLE PRECISION, INTENT(INOUT) :: W(P)
   DOUBLE PRECISION  :: TIN  ! Start Time
   DOUBLE PRECISION  :: TOUT ! End Time
   ! Optional input parameters and statistics
   INTEGER,      OPTIONAL :: ICNTRL_U(20)
   DOUBLE PRECISION,  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       OPTIONAL :: ISTATUS_U(20)
   DOUBLE PRECISION, OPTIONAL :: RSTATUS_U(20)
   INTEGER,       OPTIONAL :: IERR_U

   INTEGER,       OPTIONAL :: adICNTRL_U(20)
   DOUBLE PRECISION, OPTIONAL :: adRCNTRL_U(20)
   INTEGER,       OPTIONAL :: adISTATUS_U(20)
   DOUBLE PRECISION,  OPTIONAL :: adRSTATUS_U(20)
   INTEGER,       OPTIONAL :: adIERR_U


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

!~~~> Problem is parameter-dependent
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
       TapeSize = 2*Nd
       CALL rk_AllocateTapes(CheckptEnabled,TapeSize,N)
   ELSE
       TapeSize = ICNTRL(2)
       CALL rk_AllocateTapes(CheckptEnabled,TapeSize,N)
   END IF

   !Start clock
   CALL cpu_time(T1)
   !TLM integration
   CALL RKINT(N,Y,FEVAL,                               &
           TIN,TOUT,                                   &
           Nd, Nc,                                     &
           ATOL, RTOL,                                 &
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
   !Second order adjoint integration
   CALL  RKINT_ADJDR(N, ADY, ADJRHS, JACV,                          &
           P, W, QUADRHS, DY0DP, FEVAL,                             &
           TIN, TOUT, Nc,                                           &
           ATOL, RTOL,                                              &
           adATOL, adRTOL,                                          &
           RCNTRL, ICNTRL, RSTATUS, ISTATUS, IERR,                  &
           adRCNTRL, adICNTRL, adRSTATUS, adISTATUS, adIERR)
   !Stop clock
   CALL cpu_time(T2)
   PRINT *,'Adjoint model integration took ', (T2-T1)*1000.d0 , ' ms.'

   !check for errors that might have occurred during the adjoint model integration
   IF ((adIERR .NE. 1) .OR. (IERR .NE. 1)) THEN
       PRINT *, 'RKINT_ADJDR returned with an error!'
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

  USE PARAMETERS
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, INTENT(IN) :: T1, Y1(N)
  DOUBLE PRECISION, INTENT(OUT) :: F1(N)

  !Local variables (temporary storage)
  DOUBLE PRECISION :: tmp
  INTEGER :: i

  tmp = 0.d0

  DO i=1,np-1
     tmp = tmp + T1*params(i)*params(i+1)
  END DO

  F1(1) = (1 - Y1(2)*Y1(2))*Y1(1) - Y1(2) + tmp
  F1(2) = Y1(1)
  F1(3) = Y1(1)*Y1(1) + Y1(2)*Y1(2) + tmp*tmp

END SUBROUTINE FEVAL



SUBROUTINE TLMRHS(N,T1,Y1,DY1,TLM1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! FEVAL computes the right-hand side of the forward ODE Y' = F(T,Y)
! at (T1,Y1), i.e. F1 = F(T1,Y1).
!
! Inputs:
!  N (integer) : size of the ODE system
!  T1 (double precision) : current integration time
!  Y1 (N by 1 array, double precision) : FWD model ODE numerical solution at time T1
!  DY1 (N by 1 array, double precision) : TLM ODE numerical solution at time T1
!
! Output:
!   TLM1 (N by 1 array, double precision) : TLM ODE right-hand side at time T1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE PARAMETERS
  USE UTILS
  IMPLICIT NONE

  INTEGER :: N
  DOUBLE PRECISION :: T1
  DOUBLE PRECISION :: Y1(N), DY1(N)
  DOUBLE PRECISION :: TLM1(N)

  !Temporary variables
  DOUBLE PRECISION :: tmp
  DOUBLE PRECISION, SAVE :: df1dp(np), df3dp(np)
  INTEGER :: i

  !----------------------------------------------
  ! TANGENT LINEAR AND FUNCTION STATEMENTS
  !----------------------------------------------
  tmp = 0.d0

  DO i=1,np-1
    tmp = tmp + T1*params(i)*params(i+1)
  END DO

  !TLM model RHS = (dF/dy)*dy + (dF/dp)*dp

  !TLM1 = (dF/dy)*DY1
  TLM1(1) = -DY1(2)*(1+2*Y1(2)*Y1(1)) + DY1(1)*(1-Y1(2)*Y1(2))
  TLM1(2) = DY1(1)
  TLM1(3) = 2*DY1(2)*Y1(2) + 2*DY1(1)*Y1(1)

  !TLM1 = TLM1 + (dF/dp)*DP
  df1dp(1) = T1*params(2)

  DO i=2,np-1
     df1dp(i) = T1*params(i-1) + T1*params(i+1)
  END DO

  df1dp(np) = T1*params(np)

  df3dp = 2*tmp*df1dp

  TLM1(1) = TLM1(1) + WDOT(np,df1dp,1,DP,1)
  TLM1(3) = TLM1(3) + WDOT(np,df3dp,1,DP,1)

END SUBROUTINE TLMRHS



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
  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, INTENT(IN) :: T1, Y1(N), ADY1(N)
  DOUBLE PRECISION, INTENT(OUT) :: ADJ1(N)

  ADJ1(1) = -(1.d0 - Y1(2)*Y1(2))*ADY1(1) - ADY1(2) - 2.d0*Y1(1)*ADY1(3)
  ADJ1(2) = (1.d0 + 2.d0*Y1(1)*Y1(2))*ADY1(1) - 2*Y1(2)*ADY1(3)
  ADJ1(3) = 0.d0

END SUBROUTINE ADJRHS



SUBROUTINE JACV(N, T1, Y1, DY1, JAC1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Evaluates the RHS of the quadrature ODE system.
!  The code was generated using the TAMC automatic differentiation tool
!  (by Ralf Giering, www.autodiff.org/tamc)
!
! Inputs:
!      N (integer) : size of the forward/adjoint ODE system
!      P (integer) : size of the quadrature ODE system
!      T1 (double precision) : current integration time
!      Y1 (N by 1 array, double precision) : forward ODE numerical solution at time T1
!    DY1 (N by 1 array, double precision) : TLM ODE numerical solution at time T1
!
! Output:
!    JAC1 (N by 1 array, double precision) : Jacobian-vector product (dF/dy)*DY1 at (T1,Y1)
!                                          This is useful when performing Hermite interpolation
!                                          of the forward model trajectory.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, INTENT(IN) :: T1, Y1(N), DY1(N)
  DOUBLE PRECISION, INTENT(OUT) :: JAC1(N)

  !JAC1 = (dF/dy)*DY1
  JAC1(1) = (-(DY1(2)*(1+2*Y1(2)*Y1(1))))+DY1(1)*(1-Y1(2)*Y1(2))
  JAC1(2) = DY1(1)
  JAC1(3) = 2*DY1(2)*Y1(2)+2*DY1(1)*Y1(1)

END SUBROUTINE JACV



SUBROUTINE QUADRHS(N, P, T1, Y1, ADY1, W1)
 
   USE PARAMETERS
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: N, P
   DOUBLE PRECISION, INTENT(IN) :: T1
   DOUBLE PRECISION, INTENT(IN) :: Y1(N), ADY1(N)
   DOUBLE PRECISION, INTENT(OUT) :: W1(P)
! 
   !temporary variables
   INTEGER :: i
   DOUBLE PRECISION :: tmp
 
   tmp = 0.d0
 
   DO i=1,P-1
        tmp = tmp + T1*params(i)*params(i+1)
   END DO
! 
   ! W1 = ((dF/dp)^T * ADY1) (T1)
   W1(1) = T1*params(2)*ADY1(1) + 2.d0*T1*params(2)*ADY1(3)
   DO i=2,P-1
      W1(i) = ADY1(1) * (T1*params(i-1) + T1*params(i+1)) + &
            & 2.d0*tmp*(T1*params(i-1) + T1*params(i+1)*ADY1(3))
   END DO
   W1(P) = T1*params(P-1)*ADY1(1) + 2.d0*tmp*T1*params(P-1)*ADY1(3)
 
 END SUBROUTINE QUADRHS



SUBROUTINE DY0DP(N,P,T0,ADY0,GRAD0)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Computes the final update at time T0 to the gradient dG/dp
!
!Inputs:
!  N (integer) : size of the forward ODE system
!  P (integer) : number of parameters in the system
!  T0 (double precision) : initial time of integration
!  ADY0 (N by 1 array, double precision) : the first order adjoint solution at T0
!
!Output:
!  GRAD0 = Final correction to the Hessian-vector product (d^2G/dp^2)*(DP) at time T0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N, P
  DOUBLE PRECISION, INTENT(IN) :: T0
  DOUBLE PRECISION, INTENT(IN)  :: ADY0(N)
  DOUBLE PRECISION, INTENT(OUT) :: GRAD0(P)

  GRAD0(:) = 0.d0

END SUBROUTINE DY0DP
