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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Second order adjoint model solver for the parameter-dependent Van der Pol system
!
! Cost functional:
!                g(tF,p) = y3(tF,p)
!
! Parameters:
!       np = number of parameters
!       params(i) = 1.d0 / np
!
! Perturbation in the parameters:
!       DP(i) = 1.d0/i
!
! The code computes an approximation to the Hessian-vector product:
!
!                     W = (d^2g/dp^2) * DP
!
! using the second order adjoint method. The accuracy of the approximation can be
! controlled via the user-chosen integration tolerances, unlike the finite difference
! method where it is dependent on the finite difference parameter (denoted by "eps"
! in dr_erk_adj_chk_optc.f90).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  USE PARAMETERS
  USE UTILS

  IMPLICIT NONE

  INTEGER, PARAMETER :: N = 3
  INTEGER, PARAMETER :: P = np
  DOUBLE PRECISION :: rstatus(20), rcntrl(20)
  INTEGER :: ierr, istatus(20), icntrl(20)
  DOUBLE PRECISION :: ad_rstatus(20), ad_rcntrl(20)
  INTEGER :: ad_ierr, ad_istatus(20), ad_icntrl(20)
  DOUBLE PRECISION :: Y(N), DY(N), ADY(N), AD2Y(N)
  DOUBLE PRECISION :: W(P)
  INTEGER :: i

  !initial and final integration time
  DOUBLE PRECISION :: t0
  DOUBLE PRECISION :: tF

  INTRINSIC SQRT, REAL

  INTERFACE

    SUBROUTINE INTEGRATE(N, P, Y, DY, DP, ADY, AD2Y, W2, TIN, TOUT, &
            ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U, &
            adICNTRL_U, adRCNTRL_U, adISTATUS_U, adRSTATUS_U, adIERR_U)

   	INTEGER :: N
   	INTEGER :: P
   	DOUBLE PRECISION, INTENT(INOUT) :: Y(N), ADY(N), AD2Y(N)
   	DOUBLE PRECISION, INTENT(INOUT) :: DY(N)
   	DOUBLE PRECISION, INTENT(INOUT) :: DP(P), W2(P)
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

! Set initial time and final time
  t0 = 0.d0
  tF = 5.d0

! Set initial conditions for the forward system
  Y(1) = 0.d0
  Y(2) = 1.d0
  Y(3) = 0.d0

! Initialize second order adjoint variables
! AD2Y(tF) = (g_yy * dy + g_yp * dp)(tF) = 0.d0
  AD2Y = 0.d0

! Initialize first order adjoint variables
! ADY(t0) = g_y(tF,y,p) = (0,0,1)^T
  ADY(1) = 0.d0
  ADY(2) = 0.d0
  ADY(3) = 1.d0

! Set TLM initial condition DY = (dy0/dp)*DP = 0
! since the initial conditions of the forward model do not depend on
! the parameters
  DY = 0.d0

! Initial condition for the second order quadrature equations
! W = 0
  CALL SET2ZERO(P,W)

! Initialize the parameters
  DO i = 1,P
      params(i) = 1.d0 / P
  END DO

! Set the perturbation vector
  DO i=1,P
      DP(i) = 1.d0 / i
  END DO

! Set all integrator parameters to zero
  rcntrl(:) = 0.d0
  ad_rcntrl(:) = 0.d0
  icntrl(:) = 0
  ad_icntrl(:) = 0

!~~~> Integrate the second order adjoint system. The approximation to
! the Hessian-vector product is returned in W. One can control the accuracy
! of this approximation via the integrator tolerances (unlike in the case
! of a finite-difference approximation - see dr_erk_adj_optc.f90 - where the
! accuracy depends on the choice of eps - the finite difference coefficient)

  CALL INTEGRATE(N, P, Y, DY, DP, ADY, AD2Y, W, t0, tF, &
        & icntrl, rcntrl, istatus, rstatus, ierr, &
        & ad_icntrl, ad_rcntrl, ad_istatus, ad_rstatus, ad_ierr)

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
  PRINT *,'Second order adjoint solution:'
  PRINT *,'AD2Y(1) = ',  AD2Y(1)
  PRINT *,'AD2Y(2) = ',  AD2Y(2)
  PRINT *,'AD2Y(3) = ',  AD2Y(3)

!--- PRINT STATISTICS
  PRINT *,'No. of Jacobian-(transpose)-vector products: ', ad_istatus(1)
  PRINT *,'No. of steps (ADJ): ',              ad_istatus(2)
  PRINT *,'No. of accepted steps: ',           ad_istatus(3)
  PRINT *,'No. of rejected steps: ',           ad_istatus(4)
  PRINT *,'No. of second derivative evaluations (Hermite):', ad_istatus(5)

!   PRINT *, '****************'
!   PRINT *, 'The Hessian-vector product W = (d2g/dp2)*DP: '
!   DO i=1,P
!      PRINT *, 'W(', i, ')= ', W(i)
!   END DO

END PROGRAM MAIN



SUBROUTINE INTEGRATE(N, P, Y, DY, DP, ADY, AD2Y, W2, TIN, TOUT, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U, &
  adICNTRL_U, adRCNTRL_U, adISTATUS_U, adRSTATUS_U, adIERR_U)

   USE RKTLM_MODULE
   USE RKSOADR_MODULE
   IMPLICIT NONE
   EXTERNAL FEVAL, ADJRHS, JACV, D2Y0DP2, QUAD2RHS, SOARHS, TLMRHS

   INTEGER  :: N
   INTEGER  :: P
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N), ADY(N), AD2Y(N)
   DOUBLE PRECISION, INTENT(INOUT) :: DY(N)
   DOUBLE PRECISION, INTENT(INOUT) :: DP(P), W2(P)
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

   DOUBLE PRECISION :: ATOL(2*N), RTOL(2*N)
   DOUBLE PRECISION :: STEPMIN
   DOUBLE PRECISION :: adATOL(2*N+P), adRTOL(2*N+P)
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
       CALL rktlm_AllocateTapes(CheckptEnabled,TapeSize,N)
   ELSE
       TapeSize = ICNTRL(2)
       CALL rktlm_AllocateTapes(CheckptEnabled,TapeSize,N)
   END IF

   !Start clock
   CALL cpu_time(T1)
   !TLM integration
   CALL RKINT_TLM(N,Y,FEVAL,DY,TLMRHS,                 &
           TIN,TOUT,                                   &
           Nd, Nc,                                     &
           ATOL, RTOL,                                 &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
   !Stop clock
   CALL cpu_time(T2)
   PRINT *,'Forward + tangent linear model integration took ', (T2-T1)*1000.d0 , ' ms.'

   ! if optional parameters are given for output they 
   ! are updated with the return information
   ! update statistics for FWD integration and reset ISTATUS
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   ISTATUS(:) = 0

   PRINT *, 'Y= ', Y
   PRINT *, 'DY= ', DY

   !check for errors during the adjoint model integration
   IF (IERR .NE. 1) THEN
       PRINT *, 'RKINT_TLM returned with an error!'
       CALL rktlm_DeallocateTapes
       RETURN
   END IF

   IF (CheckptEnabled) THEN
       PRINT *, 'Wrote ', Nc, ' checkpoint(s) during FWD integration'
   END IF

   !Start clock
   CALL cpu_time(T1)
   !Second order adjoint integration
   CALL  RKINT_SOADR(N, ADY, ADJRHS, JACV,                          &
           AD2Y, SOARHS, FEVAL, TLMRHS,                             &
           P, W2, QUAD2RHS, D2Y0DP2, DP,                            &
           TIN, TOUT, Nc,                                           &
           ATOL, RTOL,                                              &
           adATOL, adRTOL,                                          &
           RCNTRL, ICNTRL, RSTATUS, ISTATUS, IERR,                  &
           adRCNTRL, adICNTRL, adRSTATUS, adISTATUS, adIERR)
   !Stop clock
   CALL cpu_time(T2)
   PRINT *,'Adjoint model integration took ', (T2-T1)*1000.d0 , ' ms.'

   !check for errors during the adjoint model integration
   IF (adIERR .NE. 1) THEN
       PRINT *, 'RKINT_SOADR returned with an error!'
       CALL rktlm_DeallocateTapes
       RETURN
   END IF

   CALL rktlm_DeallocateTapes

   STEPMIN = RSTATUS(Nhexit_tlm)
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



SUBROUTINE SOARHS(N, P, t1, y1, g_y1, g_p, ady1, g_ady1, g_adj1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Evaluates the RHS of the SOA equation
!  The code was generated using the TAMC automatic differentiation tool
!  (by Ralf Giering, www.autodiff.org/tamc)
!
! Inputs:
!      N (integer) : size of the forward/adjoint ODE system
!      P (integer) : number of system parameters
!      t1 (double precision) : current integration time
!      y1 (N by 1 array, double precision) : forward ODE numerical solution at t1
!     g_y1 (N by 1 array, double precision) : TLM ODE numerical solution at t1
!     g_p (P by 1 array, double precision) : perturbation in the parameters (delta_p)
!     ady1 (N by 1 array, double precision) : solution of the first order adjoint model at t1
!   g_ady1 (N by 1 array, double precision) : solution of the SOA model at T1
!   g_adj1 (N by 1 array, double precision) : right-hand side of the SOA equation at T1
!
! Output:
!    g_adj1 (N by 1 array, double precision) : right hand side of the SOA equation at T1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!***************************************************************
!***************************************************************
!** This routine was generated by the                         **
!** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
!***************************************************************
!***************************************************************
!==============================================
! all entries are defined explicitly
!==============================================
  IMPLICIT NONE

!==============================================
! define arguments
!==============================================
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(IN) :: P
  DOUBLE PRECISION, INTENT(IN) :: t1
  DOUBLE PRECISION, INTENT(IN) :: g_y1(N)
  DOUBLE PRECISION, INTENT(IN) :: g_p(P)
  DOUBLE PRECISION, INTENT(IN) :: y1(N)
  DOUBLE PRECISION, INTENT(IN) :: ady1(N)
  DOUBLE PRECISION, INTENT(IN) :: g_ady1(N)
  DOUBLE PRECISION, INTENT(OUT) :: g_adj1(N)


  !----------------------------------------------
  ! TANGENT LINEAR AND FUNCTION STATEMENTS
  !----------------------------------------------
  g_adj1(1) = -2*g_ady1(3)*y1(1) - g_ady1(2)           &
             - g_ady1(1)*(1.d0-y1(2)*y1(2))                &
             + 2*g_y1(2)*y1(2)*ady1(1) - 2*g_y1(1)*ady1(3)

  g_adj1(2) = -2*g_ady1(3)*y1(2) + g_ady1(1)*(1+2.d0*y1(1)*y1(2))    & 
             + g_y1(2)*(2*y1(1)*ady1(1)-2*ady1(3)) + 2*g_y1(1)*y1(2)*ady1(1)

  g_adj1(3) = 0.d0

END SUBROUTINE SOARHS



SUBROUTINE QUAD2RHS(N, P, t1, y1, g_y1, g_p, ady1, g_ady1, g_w1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Evaluates the RHS of the quadrature systems associated with the SOA
!  The code was generated using the TAMC automatic differentiation tool
!  (by Ralf Giering, www.autodiff.org/tamc)
!
! Inputs:
!      N (integer) : size of the forward/adjoint ODE system
!      P (integer) : size of the quadrature system
!      t1 (double precision) : current integration time
!      y1 (N by 1 array, double precision) : forward ODE numerical solution at t1
!     g_y1 (N by 1 array, double precision) : TLM ODE numerical solution at t1
!     g_p (P by 1 array, double precision) : perturbation in the parameters (delta_p)
!     ady1 (N by 1 array, double precision) : solution of the first order adjoint model at t1
!   g_ady1 (N by 1 array, double precision) : solution of the SOA model at T1
!
! Output:
!    g_w1 (N by 1 array, double precision) : right hand side of the quadrature equations at T1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!***************************************************************
!***************************************************************
!** This routine was generated by the                         **
!** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
!***************************************************************
!***************************************************************
  USE PARAMETERS
!==============================================
! all entries are defined explicitly
!==============================================
  IMPLICIT NONE

!==============================================
! define arguments
!==============================================
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(IN) :: P
  DOUBLE PRECISION, INTENT(IN) :: t1
  DOUBLE PRECISION, INTENT(IN) :: y1(N)
  DOUBLE PRECISION, INTENT(IN) :: g_y1(N)
  DOUBLE PRECISION, INTENT(IN) :: g_p(P)
  DOUBLE PRECISION, INTENT(IN) :: ady1(N)
  DOUBLE PRECISION, INTENT(IN) :: g_ady1(N)
  DOUBLE PRECISION, INTENT(OUT) :: g_w1(P)

  !temporary variables
  INTEGER :: i
  DOUBLE PRECISION :: tmp, g_tmp

  g_tmp = 0.d0
  tmp = 0.d0

  do i = 1, p-1
    g_tmp = g_p(i+1)*t1*params(i) + g_p(i)*t1*params(i+1) + g_tmp
    tmp = tmp + t1*params(i)*params(i+1)
  end do

  g_w1(1) = 2*g_ady1(3)*t1*params(2) + g_ady1(1)*t1*params(2) &
           + g_p(2)*(t1*ady1(1)+2*t1*ady1(3))
!  w1(1) = t1*params(2)*ady1(1) + 2.d0*t1*params(2)*ady1(3)

  do i = 2, p-1
    g_w1(i) = 2.d0*g_ady1(3)*tmp*t1*params(i+1)               &
            + g_ady1(1)*t1*(params(i-1)+params(i+1))          &
            + g_p(i-1)*(ady1(1)*t1+2.d0*tmp*t1)          &
            + g_p(i+1)*(ady1(1)*t1+2.d0*tmp*t1*ady1(3))  &
            + 2*g_tmp*t1*(params(i-1)+params(i+1)*ady1(3))

!    w1(i) = ady1(1)*(t1*params(i-1)+t1*params(i+1))+2.d0*tmp*(t1*params(i-1)+t1*params(i+1)*ady1(3))
  end do

  g_w1(p) = 2*g_ady1(3)*tmp*t1*params(p-1)                &
           + g_ady1(1)*t1*params(p-1)                     &
           +g_p(p-1)*(t1*ady1(1)+2*tmp*t1*ady1(3))   &
           +2*g_tmp*t1*params(p-1)*ady1(3)

!  w1(p) = t1*params(p-1)*ady1(1) + 2.d0*tmp*t1*params(p-1)*ady1(3)

END SUBROUTINE QUAD2RHS



SUBROUTINE D2Y0DP2(N,P,DP,T0,ADY0,AD2Y0,HESS0)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Computes the final update at time T0 to the Hessian-vector product (d^2G/dp^2)*(DP).
!
!Inputs:
!  N (integer) : size of the forward ODE system
!  P (integer) : number of parameters in the system
!  T0 (double precision) : initial time of integration
!  DP (P by 1 array, double precision) : perturbation in the parameters
!  ADY0 (N by 1 array, double precision) : the first order adjoint solution at T0
!  AD2Y0 (N by 1 array, double precision) : the second order adjoint solution at T0
!
!Output:
!  HESS0 = Final correction to the Hessian-vector product (d^2G/dp^2)*(DP) at time T0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N, P
  DOUBLE PRECISION, INTENT(IN) :: T0
  DOUBLE PRECISION, INTENT(IN) :: DP(P)
  DOUBLE PRECISION, INTENT(IN)  :: ADY0(N)
  DOUBLE PRECISION, INTENT(IN)  :: AD2Y0(N)
  DOUBLE PRECISION, INTENT(OUT) :: HESS0(P)

  HESS0(:) = 0.d0

END SUBROUTINE D2Y0DP2