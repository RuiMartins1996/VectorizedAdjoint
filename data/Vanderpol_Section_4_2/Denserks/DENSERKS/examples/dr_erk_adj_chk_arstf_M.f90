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
! Solves the Arenstorf problem and its first order adjoint ODE system
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  INTEGER, PARAMETER :: N = 4
  INTEGER, PARAMETER :: M = 4
  INTEGER, PARAMETER :: P = 0
  DOUBLE PRECISION :: t0, tF
  DOUBLE PRECISION :: rstatus(20), rcntrl(20)
  INTEGER :: ierr, istatus(20), icntrl(20)
  DOUBLE PRECISION :: ad_rstatus(20), ad_rcntrl(20)
  INTEGER :: ad_ierr, ad_istatus(20), ad_icntrl(20)
  DOUBLE PRECISION :: Y(N), ADY(N,M), W(P,M)

  INTERFACE

    SUBROUTINE INTEGRATE(N, M, P, Y, ADY, W, TIN, TOUT, &
                ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U, &
                adICNTRL_U, adRCNTRL_U, adISTATUS_U, adRSTATUS_U, adIERR_U)

      INTEGER :: N, M
      INTEGER :: P
      DOUBLE PRECISION, INTENT(INOUT) :: Y(N), ADY(N,M), W(P,M)
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

  !set the initial and final time
  t0 = 0.d0
  tF = 1.70652165601579625588917206249D0 !1/10th of the orbit period
  !tF = 17.0652165601579625588917206249D0 !full orbit period

  !set initial conditions
  Y(1)=0.994D0
  Y(2)=0.0D0
  Y(3)=0.0D0
  Y(4)=-2.00158510637908252240537862224D0

  ADY = 0.D0
  ADY(1,1) = 1.D0
  ADY(2,2) = 1.D0
  ADY(3,3) = 1.D0
  ADY(4,4) = 1.D0

  W = 0.d0

! Set all integrator parameters to zero
  rcntrl(:) = 0.d0
  ad_rcntrl(:) = 0.d0
  icntrl(:) = 0
  ad_icntrl(:) = 0

  CALL INTEGRATE(N, M, P, Y, ADY, W, t0, tF,                 &
        & icntrl, rcntrl, istatus, rstatus, ierr,            &
        & ad_icntrl, ad_rcntrl, ad_istatus, ad_rstatus, ad_ierr)

! --- PRINT FINAL SOLUTION
  PRINT *,'Final solution:'
  PRINT *,'Y(1) = ',  Y(1)
  PRINT *,'Y(2) = ',  Y(2)
  PRINT *,'Y(3) = ',  Y(3)
  PRINT *,'Y(4) = ',  Y(4)
! --- PRINT STATISTICS
  PRINT *,'No. of function calls: ', istatus(1)
  PRINT *,'No. of steps (FWD): '   , istatus(2)
  PRINT *,'No. of accepted steps: ', istatus(3)
  PRINT *,'No. of rejected steps: ', istatus(4)

  PRINT *,'*****************'

! --- PRINT ADJOINT SOLUTION
  PRINT *,'Adjoint solution:'
  PRINT *,'ADY(1,1:4) = ',  ADY(1,:)
  PRINT *,'ADY(2,1:4) = ',  ADY(2,:)
  PRINT *,'ADY(3,1:4) = ',  ADY(3,:)
  PRINT *,'ADY(4,1:4) = ',  ADY(4,:)
!--- PRINT STATISTICS
  PRINT *,'No. of Jacobian-(transpose)-vector products: ', ad_istatus(1)
  PRINT *,'No. of steps (ADJ): ',              ad_istatus(2)
  PRINT *,'No. of accepted steps: ',           ad_istatus(3)
  PRINT *,'No. of rejected steps: ',           ad_istatus(4)
  PRINT *,'No. of second derivative evaluations (Hermite):', ad_istatus(5)

END PROGRAM MAIN



SUBROUTINE INTEGRATE(N, M, P, Y, ADY, W, TIN, TOUT, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U, &
  adICNTRL_U, adRCNTRL_U, adISTATUS_U, adRSTATUS_U, adIERR_U)

   USE RK_MODULE
   USE RKADJDRM_MODULE
   IMPLICIT NONE
   EXTERNAL F, ADJRHS, QUADRHS, JACEVAL_Y, DY0DP

   INTEGER :: N, M
   INTEGER :: P
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N), ADY(N,M), W(P,M)
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
   DOUBLE PRECISION :: adATOL(N+P,M), adRTOL(N+P,M)
   DOUBLE PRECISION :: adSTEPMIN

   DOUBLE PRECISION :: T1, T2 !for timing purposes only
   LOGICAL :: CheckptEnabled

   INTEGER :: Nd, Nc, TapeSize
   INTEGER, PARAMETER :: RK5 = 4, &  !DOPRI5(4)
                         RK2 = 1, &  !Fehlberg RK2(3)
                         RK3 = 2, &  !RK3(2)
                         RK4 = 3, &  !RK4(3)
                         RK6 = 5, &  !RK6(5)
                         RK8 = 6     !DOPRI8(6)

   CheckptEnabled = .TRUE.
   Nd = 100 !write a checkpoint every Nd steps

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

   ICNTRL(2) = 100000
   adICNTRL(2) = 200000

   ICNTRL(3) = RK8
   adICNTRL(3) = RK8

   !default : dense output
   !adICNTRL(4) = 1 !3rd order Hermite
   !adICNTRL(4) = 2 !5th order Hermite

   ICNTRL(4) = 1 !DO CHECKPOINTING
   adICNTRL(6) = 1 !CHECKPOINTING WAS USED IN FWD MODE

   ICNTRL(5) = 1 !enable buffering
   adICNTRL(5) = 0 !no parameters

   ! --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
   RTOL(:)=1.0D-12
   ATOL(:)=RTOL(:)
   adRTOL(:,:)=1.0D-11
   adATOL(:,:)=adRTOL(:,:)

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
   CALL RKINT(N,Y,F,TIN,TOUT,      &
         Nd, Nc,                   &
         ATOL,RTOL,                &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
   !Stop clock
   CALL cpu_time(T2)
   PRINT *,'Forward model integration took ', (T2-T1)*1000.d0 , ' ms.'

   !check for errors during the adjoint model integration
   IF (IERR .NE. 1) THEN
       PRINT *, 'RKINT returned with an error!'
       CALL rk_DeallocateTapes
       RETURN
   END IF

   ! if optional parameters are given for output they 
   ! are updated with the return information
   ! update statistics for FWD integration and reset ISTATUS
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   ISTATUS(:) = 0

   IF (CheckptEnabled) THEN
       PRINT *, 'Wrote ', Nc, ' checkpoints during FWD integration'
   END IF

   !Start clock
   CALL cpu_time(T1)
   !Adjoint integration
   CALL RKINT_ADJDR_M(N, M, ADY, ADJRHS, JACEVAL_Y,          &
           P, W, QUADRHS, DY0DP, F,                          &
           TIN,TOUT, Nc,                                     &
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



SUBROUTINE F(N,T1,Y1,F1)

   IMPLICIT NONE
   INTEGER :: N
   DOUBLE PRECISION, INTENT(IN) :: T1, Y1(N)
   DOUBLE PRECISION, INTENT(OUT) :: F1(N)
   DOUBLE PRECISION r1, r2, mu, mu_hat

   !we're solving the Arenstorf orbit problem
   mu = 0.012277471d0 
   mu_hat = 1.d0 - mu 
   F1(1) = Y1(3)
   F1(2) = Y1(4)
   r1 = (Y1(1) + mu)**2 + Y1(2)**2
   r1 = r1*sqrt(r1)
   r2 = (Y1(1) - mu_hat)**2 + Y1(2)**2
   r2 = r2*sqrt(r2)
   F1(3) = Y1(1) + 2*Y1(4) - mu_hat*(Y1(1) + mu)/r1 - mu*(Y1(1) - mu_hat)/r2
   F1(4) = Y1(2) - 2*Y1(3) - mu_hat*Y1(2)/r1 - mu*Y1(2)/r2

END SUBROUTINE F



SUBROUTINE ADJRHS(N, T1, Y1, ADY1, ADJJAC1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !evaluates the RHS of the continuous adjoint equation at (t1,y1,ady1)
  !the continuous adjoint equation can be written as ady' = -(J^T)*ady - dg/dy 
  !where J = the Jacobian of f
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, INTENT(IN) :: T1, Y1(N), ADY1(N)
  DOUBLE PRECISION, INTENT(OUT) :: ADJJAC1(N)
  !nonconstant elements of the Jacobian
  DOUBLE PRECISION  df3dy1, df4dy1, df3dy2, df4dy2
  !parameters of the system
  DOUBLE PRECISION mu, mu_hat
  DOUBLE PRECISION N1, N2

  !we're solving the Arenstorf orbit problem
  mu = 0.012277471d0
  mu_hat = 1.d0 - mu

  N1 = (Y1(1) + mu)**2 + Y1(2)**2
  N2 = (Y1(1) - mu_hat)**2 + Y1(2)**2
  N1 = N1*N1*sqrt(N1)
  N2 = N2*N2*sqrt(N2)

  df3dy1 = 1.d0 - mu_hat*(Y1(2)**2 - 2.d0*(Y1(1) + mu)**2)/N1 - mu*(Y1(2)**2 - 2.d0*(Y1(1) - mu_hat)**2)/N2
  df3dy2 = 3.d0*mu_hat*(Y1(1) + mu)*Y1(2)/N1 + 3.d0*mu*(Y1(1) - mu_hat)*Y1(2)/N2

  df4dy1 = 3.d0*mu_hat*Y1(2)*(Y1(1) + mu)/N1 + 3.d0*mu*Y1(2)*(Y1(1) - mu_hat)/N2
  df4dy2 = 1.d0 - mu_hat*((Y1(1) + mu)**2 - 2.d0*Y1(2)**2)/N1 - mu*((Y1(1) - mu_hat)**2 - 2.d0*Y1(2)**2)/N2 

  ADJJAC1(1) = -df3dy1*ADY1(3) - df4dy1*ADY1(4)
  ADJJAC1(2) = -df3dy2*ADY1(3) - df4dy2*ADY1(4)
  ADJJAC1(3) = -ADY1(1) + 2.d0*ADY1(4)
  ADJJAC1(4) = -ADY1(2) - 2.d0*ADY1(3)

END SUBROUTINE ADJRHS



SUBROUTINE JACEVAL_Y(N, T1, Y1, F1, JAC1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !Evaluates the sensitivity equation (TLM) RHS at (t1,y1)
  !The TLM can be written as dy' = J*dy where J = df/dy = the Jacobian 
  !of f w.r.t. y in y' = f(t,y,p)
  !We evaluate the exact Jacobian (computed analytically for this problem)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, INTENT(IN) :: T1, Y1(N), F1(N)
  DOUBLE PRECISION, INTENT(OUT) :: JAC1(N)
  !nonconstant elements of the Jacobian
  DOUBLE PRECISION :: df3dy1, df4dy1, df3dy2, df4dy2
  !parameters of the system
  DOUBLE PRECISION :: mu, mu_hat
  DOUBLE PRECISION :: N1, N2

  !we're solving the Arenstorf orbit problem
  mu = 0.012277471d0 
  mu_hat = 1.d0 - mu 

  N1 = (Y1(1) + mu)**2 + Y1(2)**2
  N2 = (Y1(1) - mu_hat)**2 + Y1(2)**2
  N1 = N1*N1*sqrt(N1)
  N2 = N2*N2*sqrt(N2)

  df3dy1 = 1.d0 - mu_hat*(Y1(2)**2 - 2.d0*(Y1(1) + mu)**2)/N1 - mu*(Y1(2)**2 - 2.d0*(Y1(1) - mu_hat)**2)/N2
  df3dy2 = 3.d0*mu_hat*(Y1(1) + mu)*Y1(2)/N1 + 3.d0*mu*(Y1(1) - mu_hat)*Y1(2)/N2

  df4dy1 = 3.d0*mu_hat*Y1(2)*(Y1(1) + mu)/N1 + 3.d0*mu*Y1(2)*(Y1(1) - mu_hat)/N2
  df4dy2 = 1.d0 - mu_hat*((Y1(1) + mu)**2 - 2.d0*Y1(2)**2)/N1 - mu*((Y1(1) - mu_hat)**2 - 2.d0*Y1(2)**2)/N2 

  jac1(1) = F1(3)
  jac1(2) = F1(4)
  jac1(3) = df3dy1*F1(1) + df3dy2*F1(2) + 2.d0*F1(4)
  jac1(4) = df4dy1*F1(1) + df4dy2*F1(2) - 2.d0*F1(3)

END SUBROUTINE JACEVAL_Y



SUBROUTINE QUADRHS(N,P,T1,Y1,ADY1,W1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Evaluates the RHS of the quadrature ODE system. This is appended 
! to the ODE system if the user desires to compute the sensitivity of
! a particular cost function w.r.t. a set of parameters present 
! in the system. The quadrature ODE has the form:
! w' = (df/dp)^T * lambda + dg/dp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  IMPLICIT NONE
  INTEGER :: N
  INTEGER :: P
  DOUBLE PRECISION, INTENT(IN) :: T1
  DOUBLE PRECISION, INTENT(IN) :: Y1(N), ADY1(N)
  DOUBLE PRECISION, INTENT(OUT) :: W1(P)

  W1 = 0.d0

END SUBROUTINE QUADRHS



SUBROUTINE DY0DP(N,P,T0,ADY0,GRAD0)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 !Computes the final update to the gradient dG/dp
 !i.e. (lambda(t0)*s(t0))^T where s(t0)=dy0/dp (the sensitivity of the
 !initial solution w.r.t. the P parameters - a N x P matrix)
 !If the initial conditions are parameter-independent, then no correction
 !is needed (the subroutine must simply return a zero vector).
 !
 !Input:
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