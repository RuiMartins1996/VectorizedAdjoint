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
! Sample driver program for the Van der Pol system:
! y1' = y2
! y2' = ((1.d0 - y1**2)*y2 - y1) / eps
!
! eps = 0.01
! Initial time: t0 = 0
! Final time: tF = 4
! Initial condition: y1 = 2, y2 = 0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
! N (integer) : size of the ODE system
  INTEGER, PARAMETER :: N = 2
! Initial and final integration time moments
  DOUBLE PRECISION :: t0, tF
! Input/output parameters and error indicators
  DOUBLE PRECISION :: rstatus(20), rcntrl(20)
  INTEGER :: ierr, istatus(20), icntrl(20)
! ODE system state
  DOUBLE PRECISION :: Y(N)

  INTERFACE

	SUBROUTINE INTEGRATE(N, Y, TIN, TOUT,                &
  			     ICNTRL_U, RCNTRL_U,             &
  			     ISTATUS_U, RSTATUS_U, IERR_U )

		INTEGER :: N
   		DOUBLE PRECISION, INTENT(INOUT) :: Y(N)
   		DOUBLE PRECISION, INTENT(INOUT) :: TIN
   		DOUBLE PRECISION, INTENT(INOUT) :: TOUT
   		! Optional input parameters and statistics
   		INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   		DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   		INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   		DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   		INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U

   	END SUBROUTINE INTEGRATE


  END INTERFACE

! Set the initial and final time
  t0 = 0.d0
  tF = 4.d0

! Set initial conditions
  Y(1)=2.d0
  Y(2)=0.d0

! Set all integrator parameters to zero
  rcntrl(:) = 0.d0
  icntrl(:) = 0

! Call INTEGRATE subroutine. This is a convenience subroutine that 
! acts as a wrapper for RKINT
  CALL INTEGRATE(N, Y, t0, tF, icntrl, rcntrl, istatus, rstatus, ierr)

! --- PRINT FINAL SOLUTION (model state at tF if integration was successful)
  PRINT *,'Final solution:'
  PRINT *,'Y(1) = ',  Y(1)
  PRINT *,'Y(2) = ',  Y(2)
! --- PRINT STATISTICS
  PRINT *,'*****************'
  PRINT *,'No. of function calls: ', istatus(1)
  PRINT *,'No. of steps: '         , istatus(2)
  PRINT *,'No. of accepted steps: ', istatus(3)
  PRINT *,'No. of rejected steps: ', istatus(4)

END PROGRAM MAIN



SUBROUTINE INTEGRATE(N, Y, TIN, TOUT, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Sample wrapper for RKINT. Initializes the user-defined input 
!  parameters, calls RKINT and returns the integration results.
!
! Inputs:
!      N (integer) : size of the ODE system/forward model state
!      Y (N by 1 array, double precision) : ODE system state vector at TIN
!      TIN (double precision) : Initial integration time
!      TOUT (double precision) : Final integration time
!      ICNTRL_U (20 by 1 array, integer)
!      RCNTRL_U (20 by 1 array, double precision)
! Outputs
!      Y (N by 1 array, double precision) : ODE system state vector at TOUT
!      ISTATUS_U (20 by 1 array, integer)
!      RSTATUS_U (20 by 1 array, double precision)
!      IERR_U (integer) : error code on exit
!                    1 : integration completed successfully
!                    < 0 : an error occured. See the accompanying error 
!                        message for further details.
!
!   For more information on ICNTRL_U, RCNTRL_U, ISTATUS_U and RSTATUS_U,
!   please see the source code for RKINT in ../src/erk.f90
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Use the module containing the forward Runge-Kutta solver
   USE RK_MODULE
   IMPLICIT NONE
! Reference to the subroutine that computes the ODE right hand side
   EXTERNAL FEVAL

   INTEGER :: N
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N)
   DOUBLE PRECISION, INTENT(INOUT) :: TIN
   DOUBLE PRECISION, INTENT(INOUT) :: TOUT
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

!~~~> Maximum number of integration steps. Decrease/increase this at your own convenience
   INTEGER, PARAMETER :: MaxSteps = 100000
!~~~> File checkpointing flag
!~~~> If Checkpointing == TRUE, file checkpointing is enabled
!~~~> Else, file checkpointing is disabled (no file checkpoints will be written)
   LOGICAL :: Checkpointing
!~~~> Memory buffering flag: set it to enable buffering (TRUE) or disable it (FALSE)
   LOGICAL :: Buffering
!~~~> Nd = Number of time steps between two consecutive checkpoints
!~~~> Nc = Number of checkpoints written on disk during the forward integration stage
   INTEGER :: Nd, Nc

!~~~> Runge-Kutta method codes
   INTEGER, PARAMETER :: RK5 = 4, &  !DOPRI5(4)
                         RK2 = 1, &  !Fehlberg RK2(3)
                         RK3 = 2, &  !RK3(2)
                         RK4 = 3, &  !RK4(3)
                         RK6 = 5, &  !RK6(5)
                         RK8 = 6     !DOPRI8(6)

!~~~> Disable checkpointing and memory buffering (no adjoint integration is performed)
   Checkpointing = .FALSE.
   Buffering = .FALSE.
!~~~> Since checkpointing is disabled, the value of Nd is irrelevant
   Nd = 0

   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.d0
   ISTATUS(:) = 0
   RSTATUS(:) = 0.d0

!~~~> fine-tune the integrator:
!~~~> Set vector tolerances (specify tolerance for each component of Y)
   ICNTRL(1) = 0	! 0 - vector tolerances, 1 - scalars
!~~~> Set maximum number of steps
   ICNTRL(2) = MaxSteps
!~~~> Choose the Runge-Kutta method
   ICNTRL(3) = RK6 !RK6(5)

!~~~>Disable checkpointing and buffering
   IF (Checkpointing) THEN 
      ICNTRL(4) = 1 !Enable checkpointing if Checkpointing == TRUE
   END IF

   IF (Buffering) THEN
      ICNTRL(5) = 1 !Enable buffering
   END IF

!~~~> REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCES
   RTOL(:)=1.0D-12
   ATOL(:)=RTOL(:)

!~~~> If optional parameters are given, and if they are >0, 
!     then they overwrite default settings
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF

!~~~> Call integrator
   CALL RKINT(N,Y,FEVAL,TIN,TOUT,Nd,Nc,        &
           ATOL,RTOL,                          &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)

!~~~> If optional parameters are given for output,
!     they are updated with the return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

   IF (IERR < 0) THEN
     PRINT *,'Runge-Kutta: Unsuccessful exit at T=', TIN,' (IERR=',IERR,')'
   ENDIF

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

  DOUBLE PRECISION, PARAMETER :: eps = 0.01d0

!~~~> As a variation, the user can consider eps to be a parameter and compute
!     dy(tF)/dp or dG/dp for a given cost functional G(y)
  F1(1) = Y1(2)
  F1(2) = ((1.d0 - Y1(1)**2)*Y1(2) - Y1(1)) / eps

END SUBROUTINE FEVAL
