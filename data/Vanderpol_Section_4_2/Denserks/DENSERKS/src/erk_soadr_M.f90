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

MODULE RKSOADRM_MODULE

!~~~> use the memory buffers (tapes) implementation
  USE TLM_TAPES
!~~~> use the BLAS-like utility functions for vector computations
  USE UTILS
!~~~> use the RKINT_TLM and RKINT_SOA integrators
  USE RKTLM_MODULE
  USE RKSOA_MODULE
  IMPLICIT NONE
  PUBLIC
  SAVE

  CONTAINS

 SUBROUTINE RKINT_SOADR_M(N, M, ADYM, ADJ_RHS, JACV,                 &
           AD2YM, SOA_RHS, RHS, TLM_RHS,                             &
           P, W2M, QUAD2_RHS, D2Y0DP2, DP,                           &
           Tstart, Tend, Nc,                                         &
           FwdAbsTol, FwdRelTol,                                     &
           AdjAbsTol, AdjRelTol,                                     &
           FwdRCNTRL, FwdICNTRL, FwdRSTATUS, FwdISTATUS, FwdIERR,    &
           AdjRCNTRL, AdjICNTRL, AdjRSTATUS, AdjISTATUS, AdjIERR)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! User-callable wrapper subroutine for RKINT_ADJ (in rk_adj.f90).
! Implementation of 2-level checkpointing for second order adjoint integrations.
!
! Solves M adjoint systems, reusing the forward model integration data. Thus the
! cost for M adjoint solves is reduced by a considerable factor from that incurred
! when calling the pair RKINT/RKINT_ADJDR M times.
!
! RKINT_SOADR_M manages checkpointing by performing all the necessary
! checkpoint restores and forward trajectory recomputations.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!      N (integer)  = size of the (first/second) adjoint ODE system
!
!      M (integer)  = number of adjoint systems that need to be solved
!
!      ADYM (N by M array, double precision)
!                   = matrix of first order adjoint model final conditions (at T=Tend)
!
!      AD2YM (N by M array, double precision)
!                   = matrix of second order adjoint model final conditions (at T=Tend)
!
!      ADJ_RHS (external) = name of the subroutine computing the RHS of
!                 the adjoint ODE for a given adjoint state value. This subroutine must
!                 have the following form:
!                     SUBROUTINE ADJ_RHS(N, T1, Y1, ADY1, ADRHS1)
!                       **Inputs:**
!                       INTEGER :: N !costate size
!                       DOUBLE PRECISION :: T1 !current time moment
!                       DOUBLE PRECISION :: Y1(N) !forward model state at time T1
!                       DOUBLE PRECISION :: ADY1(N) !costate at time T1
!                       **Output:**
!                       DOUBLE PRECISION :: ADRHS1(N) !RHS of the adjoint ODE
!
!      SOA_RHS (external) = name of the subroutine computing the RHS of the second order adjoint
!                 ODE. The subroutine must have the following form:
!                    SUBROUTINE SOA_RHS(N,T1,Y1,DY1,DP,ADY1,AD2Y1,AD2RHS1)
!                       **Inputs:**
!                       INTEGER :: N,P !costate size and number of parameters (if any)
!                       DOUBLE PRECISION :: T1 !current time moment
!                       DOUBLE PRECISION :: Y1(N) !forward model state at T1
!                       DOUBLE PRECISION :: DY1(N) !TLM solution at T1 (y_p*dp(T1) -> see SOA equation)
!                       DOUBLE PRECISION :: DP(P) !perturbation in the parameters
!                       DOUBLE PRECISION :: ADY1(N) !first orde adjoint model state at T1
!                       DOUBLE PRECISION :: AD2Y1(N) !second order adjoint model state at T1
!                       **Output:**
!                       DOUBLE PRECISION :: AD2RHS1(N) !RHS of the SOA ODE
!
!      JACV (external) = name of subroutine computing the value of the Jacobian of
!            F times a given user vector v, i.e. Jac*v. This is needed for the 5th
!            order Hermite interpolation in the adjoint computation. In this case,
!            the procedure is used to compute an approximation of the second
!            derivative of the forward solution: y''(t). If 3rd order
!            Hermite interpolation or dense output are used, the user can pass a
!            dummy pointer for this argument. Otherwise, the subroutine must have the
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
!      P  (integer) = size of the quadrature ODE system (i.e. number of parameters)
!
!      W2M (P by M array, double precision) = quadrature ODE system states at (T->Tend)
!
!      QUAD2_RHS (external) = name of the subroutine computing the RHS of the quadrature ODE
!                  associated with the second-order adjoint system. If the user wants to
!                  compute information about the Hessian of G w.r.t. the initial conditions,
!                  then he or she can pass in a dummy pointer for this parameter (quadrature
!                  equations are not present in the system). Otherwise, the subroutine must have
!                  the following form:
!                        QUAD2_RHS(N,P,T1,Yval,DYVal,DP,ADY1,AD2Y1,Quad)
!                       **Inputs:**
!                       INTEGER :: N,P !costate size and number of parameters (if any)
!                       DOUBLE PRECISION :: T1 !current time moment
!                       DOUBLE PRECISION :: Yval(N) !forward model state at T1
!                       DOUBLE PRECISION :: DYval(N) !TLM model state at T1
!                       DOUBLE PRECISION :: DP(P) !perturbation in the parameters
!                       DOUBLE PRECISION :: ADY1(N) !first order adjoint model state at T1
!                       DOUBLE PRECISION :: AD2Y1(N) !second order forward model state at T1
!                       **Output:**
!                       DOUBLE PRECISION :: Quad(P) !RHS of the quadrature ODE
!
!      Tstart (double precision), Tend (double precision) :
!        [Tstart,Tend]  = time range of integration; RKINT_SOADR requires Tstart<=Tend
!
!      Nc (integer) = total number of checkpoints that have written during the forward
!                 model integration.
!                (Nc = 0 if checkpointing was disabled)
!
!      D2YODP2 (external) = name of the subroutine computing the final correction to
!                  the Hessian-vector product Hess*DP at time T->Tstart. This is needed if
!                  the initial conditions of the tangent linear model depend on the system
!                  parameters. In this case, it is necessary to add the term (y_p^T)*AD2Y(Tstart)
!                  to W2 at T->Tstart. If no parameters are present, the user can
!                  pass in a dummy pointer. Otherwise, the following subroutine header is
!                  required:
!                     SUBROUTINE D2Y0DP2(N,P,T0,DP,ADY0,AD2Y0,W20)
!                       **Inputs:**
!                       INTEGER :: N,P !costate size and number of parameters (if any)
!                       DOUBLE PRECISION :: T0       !initial integration time
!                       DOUBLE PRECISION :: DP(P)    !perturbation in the parameters
!                       DOUBLE PRECISION :: ADY0(N)  !adjoint model state at Tstart
!                       DOUBLE PRECISION :: AD2Y0(N) !adjoint model state at Tstart
!                       **Output:**
!                       DOUBLE PRECISION :: W20(P) !Hessian-vector product correction
!
!      RHS (external) = name of subroutine computing the value of F(T,Y)
!                  SUBROUTINE RHS(N,T,Y,F)
!                  DOUBLE PRECISION T, Y(N),F(N)
!                  F(1) = .... etc.
!
!     FwdRelTol (2*N by M array, double precision),
!     FwdAbsTol (2*N by M array, double precision) :
!                                  = user precribed relative and absolute
!                                  accuracy for the forward model and TLM integration; these are
!                                  needed in case forward/TLM trajectory recomputations need to
!                                  be performed.
!                       FwdRelTol(1:N),FwdAbsTol(1:N) = relative/absolute tolerances used when
!                                     the forward model was integrated; needed in case forward trajectory
!                                     recomputations are necessary
!                       FwdRelTol(N+1:2*N),FwdAbsTol(N+1:2*N) = relative/absolute tolerances used
!                                     when the TLM model was integrated; needed in case TLM trajectory
!                                     recomputations are necessary
!
!     AdjRelTol (2*N+P by M array, double precision),
!     AdjAbsTol (2*N+P by M array, double precision)
!                                    = user prescribed relative and absolute
!                                      accuracy for the adjoint model integration
!                       AdjRelTol(1:N,i),AdjAbsTol(1:N,i) = relative/absolute tolerances for
!                                     the first order adjoint model variables ADYM(:,i)
!                       AdjRelTol(N+1:2*N,i),AdjAbsTol(N+1:2*N,i) = relative/absolute tolerances for
!                                     the second order adjoint model variables AD2YM(:,i)
!                       AdjRelTol(2*N+1:P,i),AdjAbsTol(2*N+1:P,i) = relative/absolute tolerances for
!                                     the quadrature variables W2M(:,i)
!
!     Fwd/AdjICNTRL (20 by 1 array, integer) = integer input parameters
!                                             (see rk_tlm.f90 and rk_soa.f90)
!     Fwd/AdjRCNTRL (20 by 1 array, double precision) = real input parameters
!                                             (see rk_tlm.f90 and rk_soa.f90)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     OUTPUT ARGUMENTS:
!
!     ADYM (N by M array, double precision)  -> matrix of first order adjoint model states (at T->Tstart)
!     AD2YM (N by M array, double precision) -> matrix of second order adjoint model states (at T->Tstart)
!     W2M (P by M array, double precision)   -> Hessian-vector product (d2G/dp2)*DP
!
!-    Fwd/AdjISTATUS(1:20)   -> integer output parameters
!-    Fwd/AdjRSTATUS(1:20)   -> real output parameters
!-    Fwd/AdjIERR            -> job status upon return
!                         success (positive value) or
!                         failure (negative value)
!
!   (See rk.f90, rk_adj.f90, rk_soa.f90 for more details on these parameters)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   IMPLICIT NONE
   EXTERNAL ADJ_RHS, SOA_RHS, JACV
   EXTERNAL QUAD2_RHS, RHS, TLM_RHS
   EXTERNAL D2Y0DP2
   INTEGER,          INTENT(INOUT) :: N
   INTEGER,          INTENT(INOUT) :: Nc
   INTEGER,          INTENT(INOUT) :: M
   INTEGER,          INTENT(INOUT) :: P
   DOUBLE PRECISION, INTENT(INOUT) :: ADYM(N,M), AD2YM(N,M)
   DOUBLE PRECISION, INTENT(INOUT) :: W2M(P,M), DP(P)
   DOUBLE PRECISION, INTENT(INOUT) :: Tstart,Tend
   DOUBLE PRECISION, INTENT(INOUT) :: FwdAbsTol(2*N), FwdRelTol(2*N)
   DOUBLE PRECISION, INTENT(INOUT) :: AdjAbsTol(2*N+P,M), AdjRelTol(2*N+P,M)
   INTEGER,          INTENT(INOUT) :: FwdICNTRL(20), AdjICNTRL(20)
   DOUBLE PRECISION, INTENT(INOUT) :: FwdRCNTRL(20), AdjRCNTRL(20)
   INTEGER,          INTENT(INOUT) :: FwdISTATUS(20), AdjISTATUS(20)
   DOUBLE PRECISION, INTENT(INOUT) :: FwdRSTATUS(20), AdjRSTATUS(20)
   INTEGER,          INTENT(OUT)   :: FwdIERR, AdjIERR

   DOUBLE PRECISION, PARAMETER :: ZERO = 0.d0
   DOUBLE PRECISION :: chk_H, chk_T, chk_T_Next
   DOUBLE PRECISION :: chk_Y(N), chk_DY(N), W2update(P)
! FWD mode statistics:
   DOUBLE PRECISION :: FwdISTATUS_T(20)
! ADJ mode statistics:
   DOUBLE PRECISION :: AdjISTATUS_T(20)
   INTEGER :: Nd = 0
   INTEGER :: idx, i
   LOGICAL :: Checkpointing

!~~~> AdjICNTRL(6) = 1 if we employ checkpointing
!                  = 0 otherwise
   IF (AdjICNTRL(6) == 0) THEN
      Checkpointing = .FALSE.
   ELSEIF (AdjICNTRL(6) > 0) THEN
      Checkpointing = .TRUE.
   ELSE
      PRINT * ,'Invalid checkpointing option: AdjICNTRL(6)=', AdjICNTRL(6)
      CALL rkdr_ErrorMsg(-1,Tstart,ZERO,AdjIERR)
      RETURN
   END IF

   FwdISTATUS_T(:) = 0
   AdjISTATUS_T(:) = 0

   IF (.NOT. Checkpointing) THEN

       DO i=1,M

          !all the necessary data for interpolation is kept in the memory tapes
          !just call RKINT_SOA to integrate the adjoint system from Tend to Tstart
          CALL RKINT_SOA(N, ADYM(:,i), ADJ_RHS, JACV,                   &
                         AD2YM(:,i), SOA_RHS,                            &
                         P, W2M(:,i), QUAD2_RHS, DP,                    &
                         Tstart, Tend,                                  &
                         AdjAbsTol(:,i), AdjRelTol(:,i),                &
                         AdjRCNTRL, AdjICNTRL, AdjRSTATUS, AdjISTATUS, AdjIERR)

          !reinitialize memory buffer pointers to prepare for next adjoint integration
          CALL tlm_tapes_ReinitPointers

          AdjISTATUS_T(1:6) = AdjISTATUS_T(1:6) + AdjISTATUS(1:6)

          IF (AdjIERR .LE. 0) THEN
                !an error occured, no point to continue the integration
                RETURN
          END IF

       END DO

   ELSE

      !we chose to checkpoint data during the FWD integration
      !how many checkpoints were written? (at least one)
      IF (Nc .EQ. 1) THEN

          DO i=1,M

               CALL RKINT_SOA(N, ADYM(:,i), ADJ_RHS, JACV,                &
                          AD2YM(:,i), SOA_RHS,                            &
                          P, W2M(:,i), QUAD2_RHS, DP,                     &
                          Tstart, Tend,                                   &
                          AdjAbsTol(:,i), AdjRelTol(:,i),                 &
                          AdjRCNTRL, AdjICNTRL, AdjRSTATUS, AdjISTATUS, AdjIERR)

              AdjISTATUS_T(1:6) = AdjISTATUS_T(1:6) + AdjISTATUS(1:6)

              !reinitialize memory buffer pointers to prepare for next adjoint integration
              CALL tlm_tapes_ReinitPointers

              IF (AdjIERR .LE. 0) THEN
                  !an error occured, no point to continue the integration
                  RETURN
              END IF

          END DO

      ELSE !(Nc > 1)

          !read checkpoints in reverse order
          CALL tlm_tapes_ReadCheckpt(Nc,chk_T,chk_H,chk_Y,chk_DY)

          DO i=1,M

              !the first time there's no need to recreate FWD trajectory for [chk_T,Tend]
               CALL RKINT_SOA(N, ADYM(:,i), ADJ_RHS, JACV,                &
                          AD2YM(:,i), SOA_RHS,                            &
                          P, W2M(:,i), QUAD2_RHS, DP,                     &
                          chk_T, Tend,                                    &
                          AdjAbsTol(:,i), AdjRelTol(:,i),                 &
                          AdjRCNTRL, AdjICNTRL, AdjRSTATUS, AdjISTATUS, AdjIERR)

               !reinitialize memory buffer pointers to prepare for next adjoint integration
               CALL tlm_tapes_ReinitPointers

               AdjISTATUS_T(1:6) = AdjISTATUS_T(1:6) + AdjISTATUS(1:6)

               IF (AdjIERR .LE. 0) THEN
                   !an error occured, no point to continue the integration
                   RETURN
               END IF

          END DO

          chk_T_Next = chk_T

          !read the rest of the checkpoints
          !for these subintervals, the fwd trajectory needs to be recalculated
          DO idx = Nc-1,1,-1

             CALL tlm_tapes_ReadCheckpt(idx,chk_T,chk_H,chk_Y,chk_DY)
             CALL tlm_tapes_Flush

             !set the starting time step to allow a hot restart
             FwdRCNTRL(3) = chk_H
             FwdICNTRL(4) = 0 !not checkpointing anymore
             FwdICNTRL(5) = 1 !but we still need the memory buffers

             CALL RKINT_TLM(N, chk_Y, RHS, chk_DY, TLM_RHS,           &
                    chk_T, chk_T_Next,                                &
                    Nd, Nc,                                           &
                    FwdAbsTol, FwdRelTol,                             &
                    FwdRCNTRL,FwdICNTRL,FwdRSTATUS,FwdISTATUS,FwdIERR)

             FwdISTATUS_T(1:4) = FwdISTATUS_T(1:4) + FwdISTATUS(1:4)

             IF (FwdIERR .LE. 0) THEN
                  STOP 1
             END IF

             DO i=1,M

                  CALL RKINT_SOA(N, ADYM(:,i), ADJ_RHS, JACV,        &
                            AD2YM(:,i), SOA_RHS,                     &
                            P, W2M(:,i), QUAD2_RHS, DP,              &
                            chk_T, chk_T_Next,                       &
                            AdjAbsTol(:,i), AdjRelTol(:,i),          &
                            AdjRCNTRL, AdjICNTRL, AdjRSTATUS, AdjISTATUS, AdjIERR)

                 !reinitialize memory buffer pointers to prepare for next adjoint integration
                 CALL tlm_tapes_ReinitPointers

                 AdjISTATUS_T(1:6) = AdjISTATUS_T(1:6) + AdjISTATUS(1:6)

                 IF (AdjIERR .LE. 0) THEN
                      RETURN
                 END IF

             END DO

             chk_T_Next = chk_T

          END DO

      END IF

   END IF

!~~~> Hessian-vector update at T0
   DO i=1,M

      IF (AdjICNTRL(5) .EQ. 1) THEN

          CALL D2Y0DP2(N,P,DP,Tstart,ADYM(:,i),AD2YM(:,i),W2update)
          CALL WAXPY(P,1.d0,W2update,1,W2M(:,i),1)

      END IF

   END DO

   IF (Checkpointing) THEN 
      CALL tlm_tapes_Cleanup
   END IF

   !Final statistics
   AdjISTATUS(1:6) = AdjISTATUS_T(1:6)
   FwdISTATUS(1:4) = FwdISTATUS_T(1:4)

   AdjIERR = 1 !Successful integration
   FwdIERR = 1 !Successful integration

  CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE rkdr_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR

   IERR = Code
   PRINT * , &
     'Forced exit from RKINT_SOADR_M due to the following error:' 

   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper checkpointing option AdjICNTRL(6)=', AdjICNTRL(6)
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT

   PRINT *, "T=", T, "and H=", H

 END SUBROUTINE rkdr_ErrorMsg


 END SUBROUTINE RKINT_SOADR_M


END MODULE RKSOADRM_MODULE