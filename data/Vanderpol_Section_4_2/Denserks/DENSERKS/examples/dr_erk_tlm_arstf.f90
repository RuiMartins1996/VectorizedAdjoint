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
! Generic driver for the RKINT and RKINT_ADJ integrators
! By default, it solves the Arenstorf problem and its adjoint ODE system
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  INTEGER, PARAMETER :: N = 4
  DOUBLE PRECISION :: t0, tF
  DOUBLE PRECISION :: rstatus(20), rcntrl(20)
  INTEGER :: ierr, istatus(20), icntrl(20)
  DOUBLE PRECISION :: Y(N), DY(N)

  INTERFACE

	SUBROUTINE INTEGRATE(N, Y, DY, TIN, TOUT,           &
  			    ICNTRL_U, RCNTRL_U, ISTATUS_U,  &
  			    RSTATUS_U, IERR_U )

   		INTEGER :: N
   		DOUBLE PRECISION, INTENT(INOUT) :: Y(N), DY(N)
   		DOUBLE PRECISION, INTENT(INOUT) :: TIN  ! Start Time
   		DOUBLE PRECISION, INTENT(INOUT) :: TOUT ! End Time
   		! Optional input parameters and statistics
   		INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   		DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   		INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   		DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   		INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U

	END SUBROUTINE INTEGRATE

  END INTERFACE

! set the initial and final time
  t0 = 0.d0
! tF = 17.0652165601579625588917206249D0 !full orbit period
  tF = 1.70652165601579625588917206249D0 !1/10th of the orbit's period

  !set initial conditions
  Y(1)=0.994D0
  Y(2)=0.0D0
  Y(3)=0.0D0
  Y(4)=-2.00158510637908252240537862224D0

  DY(1:4) = 0.d0
  DY(1) = 1.d0

! Set all integrator parameters to zero
  rcntrl(:) = 0.d0
  icntrl(:) = 0

  CALL INTEGRATE(N, Y, DY, t0, tF, icntrl, rcntrl, istatus, rstatus, ierr)

! --- PRINT FINAL SOLUTION
  PRINT *,'Final solution:'
  PRINT *,'Y(1) = ',  Y(1)
  PRINT *,'Y(2) = ',  Y(2)
  PRINT *,'Y(3) = ',  Y(3)
  PRINT *,'Y(4) = ',  Y(4)
! --- PRINT FINAL SOLUTION
  PRINT *,'*****************'
  PRINT *,'Final TLM solution:'
  PRINT *,'DY(1) = ',  DY(1)
  PRINT *,'DY(2) = ',  DY(2)
  PRINT *,'DY(3) = ',  DY(3)
  PRINT *,'DY(4) = ',  DY(4)
! --- PRINT STATISTICS
  PRINT *,'*****************'
  PRINT *,'No. of function calls: ', istatus(1)
  PRINT *,'No. of steps: '         , istatus(2)
  PRINT *,'No. of accepted steps: ', istatus(3)
  PRINT *,'No. of rejected steps: ', istatus(4)

END PROGRAM MAIN



SUBROUTINE INTEGRATE(N, Y, DY, TIN, TOUT, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Simple wrapper for RKINT_TLM.
!  Initializes the user-defined integrator parameters,
!  calls RKINT and then RKINT_ADJDR and returns the integration results.
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
!      Y (N by 1 array, double precision) : FWD ODE system state vector at TOUT
!      DY (N by 1 array, double precision) : TLM ODE system state vector at TOUT
!      ISTATUS_U (20 by 1 array, integer)
!      RSTATUS_U (20 by 1 array, double precision)
!      IERR_U (integer) : forward model integrator error code on exit
!                    1 : integration completed successfully
!                    < 0 : an error occured. See the accompanying error
!                        message for further details.
!
!   For more information on ICNTRL_U, RCNTRL_U, ISTATUS_U and RSTATUS_U,
!   please see the source code for RKINT_TLM in ../src/erk_tlm.f90
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   USE RKTLM_MODULE
   IMPLICIT NONE
   EXTERNAL FEVAL, TLM_RHS

   INTEGER :: N
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N), DY(N)
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

   DOUBLE PRECISION :: ATOL(2*N), RTOL(2*N)
   DOUBLE PRECISION :: STEPMIN

   INTEGER :: MaxSteps = 50000
   LOGICAL :: Checkpointing
   INTEGER :: Nd, Nc

   INTEGER, PARAMETER :: RK5 = 4, &  !DOPRI5(4)
                         RK2 = 1, &  !Fehlberg RK2(3)
                         RK3 = 2, &  !RK3(2)
                         RK4 = 3, &  !RK4(3)
                         RK6 = 5, &  !RK6(5)
                         RK8 = 6     !DOPRI8(6)

   !Forward sensitivity analysis only, so no checkpointing needed
   Nd = 0
   Checkpointing = .FALSE.

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.d0
   ISTATUS(:) = 0
   RSTATUS(:) = 0.d0

   !~~~> Fine-tune the integrator:
   ICNTRL(1) = 0        ! 0 - vector tolerances, 1 - scalars
   ICNTRL(2) = MaxSteps
   ICNTRL(3) = RK8      !RK method of choice

   !~~~> Forward sensitivity analysis only, so disable memory buffering
   !and file checkpointing
   ICNTRL(4) = 0 !No checkpointing
   ICNTRL(5) = 0 !Disable buffering

   !~~~> REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCES
   RTOL(:)=1.0D-9
   ATOL(:)=RTOL(:)

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF

   CALL RKINT_TLM(N,Y,FEVAL,DY,TLM_RHS,                &
           TIN,TOUT,                                   &
           Nd, Nc,                                     &
           ATOL, RTOL,                                 &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)

   STEPMIN = RSTATUS(Nhexit_tlm)

   ! if optional parameters are given for output they
   ! are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

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
  INTEGER :: N
  DOUBLE PRECISION, INTENT(IN) :: T1, Y1(N)
  DOUBLE PRECISION, INTENT(OUT) :: F1(N)
  DOUBLE PRECISION r1, r2, mu, mu_hat

  !we're solving the Arenstorf orbit problem
  mu = 0.012277471d0 
  mu_hat = 1.d0 - mu 
  f1(1) = y1(3)
  f1(2) = y1(4)
  r1 = (y1(1) + mu)**2 + y1(2)**2
  r1 = r1*sqrt(r1)
  r2 = (y1(1) - mu_hat)**2 + y1(2)**2
  r2 = r2*sqrt(r2)
  f1(3) = y1(1) + 2*y1(4) - mu_hat*(y1(1) + mu)/r1 - mu*(y1(1) - mu_hat)/r2
  f1(4) = y1(2) - 2*y1(3) - mu_hat*y1(2)/r1 - mu*y1(2)/r2

END SUBROUTINE FEVAL


SUBROUTINE TLM_RHS(N, T1, Y1, G_Y1, G_F1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  integer, intent(in) :: N
  double precision, intent(out) :: g_f1(n)
  double precision, intent(in) :: g_y1(n)
  double precision, intent(in) :: t1
  double precision, intent(in) :: y1(n)

!==============================================
! define local variables
!==============================================
  double precision :: g_r1
  double precision :: g_r2
  double precision :: mu
  double precision :: mu_hat
  double precision :: r1
  double precision :: r2

!----------------------------------------------
! The tangent linear model was generated using TAMC
! TANGENT LINEAR MODEL STATEMENTS
!----------------------------------------------
  mu = 0.012277471d0
  mu_hat = 1.d0-mu
  g_f1(1) = g_y1(3)
  g_f1(2) = g_y1(4)
  g_r1 = 2*g_y1(2)*y1(2)+2*g_y1(1)*(y1(1)+mu)
  r1 = (y1(1)+mu)**2+y1(2)**2
  g_r1 = g_r1*(r1*(1./(2.*sqrt(r1)))+sqrt(r1))
  r1 = r1*sqrt(r1)
  g_r2 = 2*g_y1(2)*y1(2)+2*g_y1(1)*(y1(1)-mu_hat)
  r2 = (y1(1)-mu_hat)**2+y1(2)**2
  g_r2 = g_r2*(r2*(1./(2.*sqrt(r2)))+sqrt(r2))
  r2 = r2*sqrt(r2)
  g_f1(3) = g_r1*(mu_hat*(y1(1)+mu)/(r1*r1))+g_r2*(mu*(y1(1)-mu_hat)/(r2*r2))+2*g_y1(4)+g_y1(1)*(1-mu_hat/r1-mu/r2)
  g_f1(4) = g_r1*(mu_hat*y1(2)/(r1*r1))+g_r2*(mu*y1(2)/(r2*r2))-2*g_y1(3)+g_y1(2)*(1-mu_hat/r1-mu/r2)

END SUBROUTINE TLM_RHS