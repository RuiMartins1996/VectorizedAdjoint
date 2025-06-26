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

MODULE RKADJDR_MODULE

!~~~> use the memory buffers (tapes) implementation
  USE TAPES
!~~~> use the BLAS-like utility functions for vector computations
  USE UTILS
!~~~> Use the RKINT and RKINT_ADJ subroutines for model integrations
  USE RK_MODULE
  USE RKADJ_MODULE
  IMPLICIT NONE
  PUBLIC
  SAVE

 CONTAINS

 SUBROUTINE RKINT_ADJDR(N, ADY, ADJ_RHS, JACV,                      &
           P, W, QUAD_RHS, DY0DP, RHS,                              &
           Tstart, Tend, Nc,                                        &
           FwdAbsTol, FwdRelTol,                                    &
           AdjAbsTol, AdjRelTol,                                    &
           FwdRCNTRL, FwdICNTRL, FwdRSTATUS, FwdISTATUS, FwdIERR,   &
           AdjRCNTRL, AdjICNTRL, AdjRSTATUS, AdjISTATUS, AdjIERR)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! User-callable wrapper subroutine for RKINT_ADJ (in rk_adj.f90).
! Implementation of 2-level checkpointing for first order adjoint integrations.
!
! RKINT_ADJDR manages checkpointing by performing all the necessary
! checkpoint restores and forward trajectory recomputations.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!      N (integer)       = size of the adjoint ODE system
!
!      ADY (N by 1 array, double precision)  = vector of adjoint final conditions (at T=Tend)
!
!      ADJ_RHS (external) = name of the subroutine computing the RHS of
!                the adjoint ODE for a given adjoint state value. This subroutine must
!                have the following form:
!                     SUBROUTINE ADJ_RHS(N, T1, Y1, ADY1, ADJ1)
!                       **Inputs:**
!                       INTEGER :: N !costate size
!                       DOUBLE PRECISION :: T1 !current time moment
!                       DOUBLE PRECISION :: Y1(N) !forward model state at time T1
!                       DOUBLE PRECISION :: ADY1(N) !costate at time T1
!                       **Output:**
!                       DOUBLE PRECISION :: ADJ1(N) !RHS of the adjoint ODE
!
!                It is recommended that the reverse mode of automatic differentiation be
!                used to generate the code that computes the product Jac^T * lambda. This
!                approach is very efficient (avoids computing and storing the entire Jacobian
!                matrix), and, unlike a finite difference approximation, it does not incur any
!                truncation errors (the result will be exact up to working accuracy).
!
!      JACV  (external) = name of the subroutine computing the value of the Jacobian of
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
!             It is recommended that the forward mode of automatic differentiation be used
!             to generate the code for this method.
!
!      QUAD_RHS (external) = name of the subroutine computing the RHS of the quadrature
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
!                The reverse mode of automatic differentiation can (and should) be used
!                to generate the (Jacobian-transpose)-vector product required in the 
!                right hand side of the quadrature ODEs.
!
!      P  (integer)  = size of the quadrature ODE system (i.e. number of parameters)
!
!      W (P by 1 array, double precision) = quadrature ODE system state at (T->Tend)
!
!      Tstart (double precision), Tend (double precision) :
!          [Tstart,Tend]  = time range of integration. RKINT_ADJDR requires that Tstart<=Tend
!
!      Nc (integer) = total number of checkpoints that have written during the forward
!                 model integration.
!                (Nc = 0 if checkpointing was disabled)
!
!      DYODP (external) = name of the subroutine computing the final correction to
!                 the gradient dG/dp at time T->Tstart. This is needed if the initial
!                 conditions for the forward model depend on the system parameters. In
!                 this case, it is necessary to add the term (dY(t0)/dp)*ADY(Tstart) to
!                 the W at T->Tstart. If no parameters are present, the user can, again,
!                 pass a dummy pointer. Otherwise, the following subroutine header is
!                 required:
!                     SUBROUTINE DY0DP(N,P,T0,ADY0,GRAD0)
!                       **Inputs:**
!                       INTEGER :: N,P !costate size and number of parameters (if any)
!                       DOUBLE PRECISION :: T0  !initial integration time
!                       DOUBLE PRECISION :: ADY0(N) !adjoint model state at Tstart
!                       **Output:**
!                       DOUBLE PRECISION :: GRAD0(P) !gradient correction
!
!      RHS (external) = name of subroutine computing the value of F(T,Y)
!                  SUBROUTINE RHS(N,T,Y,F)
!                  DOUBLE PRECISION T, Y(N),F(N)
!                  F(1) = .... etc.
!
!     FwdRelTol (N by 1 array, double precision),
!     FwdAbsTol (N by 1 array, double precision) :
!                                  = user precribed relative and absolute
!                                  accuracy for the forward model integration; these are
!                                  needed in case forward trajectory recomputations need to
!                                  be performed
!
!     AdjRelTol (N+P by 1 array, double precision),
!     AdjAbsTol (N+P by 1 array, double precision)
!                                    = user prescribed relative and absolute
!                                      accuracy for the adjoint model integration
!                       AdjRelTol(1:N),AdjAbsTol(1:N) = relative/absolute tolerances for
!                                     the first order adjoint model variables (ADY)
!                       AdjRelTol(N+1:P),AdjAbsTol(N+1:P) = relative/absolute tolerances for
!                                     the quadrature variables (W)
!
!     Fwd/AdjICNTRL (20 by 1 array, integer) = integer input parameters (see rk.f90 and rk_adj.f90)
!     Fwd/AdjRCNTRL (20 by 1 array, double precision) = real input parameters (see rk.f90 and rk_adj.f90)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     OUTPUT ARGUMENTS:
!
!     ADY (N by 1 array, double precision) -> vector of final adjoint states (at T->Tstart), i.e. solution
!                  of the first order adjoint final value problem
!
!     W (P by 1 array, double precision) ->  solution of the quadrature ODE system with final correction (if applicable)
!                            i.e. the gradient of the cost function dg/dp
!
!-    Fwd/AdjISTATUS (20 by 1 array, integer)
!                       -> integer output parameters
!-    Fwd/AdjRSTATUS (20 by 1 array, double precision)
!                       -> real output parameters
!-    Fwd/AdjIERR (integer)
!                       -> job status upon return:
!                            success (positive value) or
!                            failure (negative value, value equals error code)
!     (See rk.f90 and rk_adj.f90 for more details on these parameters)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   IMPLICIT NONE
   EXTERNAL ADJ_RHS, JACV, QUAD_RHS, DY0DP, RHS
   INTEGER,          INTENT(IN) :: N, Nc
   INTEGER,          INTENT(IN) :: P
   DOUBLE PRECISION, INTENT(INOUT) :: ADY(N)
   DOUBLE PRECISION, INTENT(INOUT) :: W(P)
   DOUBLE PRECISION, INTENT(INOUT) :: Tstart,Tend
   DOUBLE PRECISION, INTENT(INOUT) :: FwdAbsTol(N), FwdRelTol(N)
   DOUBLE PRECISION, INTENT(INOUT) :: AdjAbsTol(N+P), AdjRelTol(N+P)
   INTEGER,          INTENT(INOUT) :: FwdICNTRL(20), AdjICNTRL(20)
   DOUBLE PRECISION, INTENT(INOUT) :: FwdRCNTRL(20), AdjRCNTRL(20)
   INTEGER,          INTENT(INOUT) :: FwdISTATUS(20), AdjISTATUS(20)
   DOUBLE PRECISION, INTENT(INOUT) :: FwdRSTATUS(20), AdjRSTATUS(20)
   INTEGER,          INTENT(OUT)   :: FwdIERR, AdjIERR

   DOUBLE PRECISION, PARAMETER :: ZERO = 0.d0
   DOUBLE PRECISION :: chk_H, chk_T, chk_T_Next
   DOUBLE PRECISION :: chk_Y(N), Wupdate(P)
! FWD:
   DOUBLE PRECISION :: FwdISTATUS_T(20)
! ADJ:
   DOUBLE PRECISION :: AdjISTATUS_T(20)
   INTEGER :: Nd = 0, Nc1 = 0
   INTEGER :: idx
   LOGICAL :: Checkpointing

   !ICNTRL(6) = 1 if we employ checkpointing
   !          = 0 otherwise
   IF (AdjICNTRL(6) == 0) THEN
      Checkpointing = .FALSE.
   ELSEIF (AdjICNTRL(6) > 0) THEN
      Checkpointing = .TRUE.
   ELSE
      PRINT * ,'Invalid checkpointing option: AdjICNTRL(6)=', AdjICNTRL(6)
      CALL rkdr_ErrorMsg(-1,Tstart,ZERO,AdjIERR)
      RETURN
   END IF

   !Check validity of Nc (must be positive if checkpointing is enabled)
   IF (Checkpointing .AND. (Nc .LE. 0)) THEN
      PRINT * ,'Invalid value Nc = ', Nc
      CALL rkdr_ErrorMsg(-2,Tstart,ZERO,AdjIERR)
      RETURN
   END IF

   FwdISTATUS_T(:) = 0
   AdjISTATUS_T(:) = 0

   IF (.NOT. Checkpointing) THEN

       !all the necessary data for interpolation is kept in the tapes (in memory)
       !just call RK_ADJ to integrate the adjoint system from Tend to Tstart
       CALL RKINT_ADJ(N, ADY, ADJ_RHS, JACV,                        &
                    P, W, QUAD_RHS,                                 &
                    Tstart, Tend,                                   &
                    AdjAbsTol, AdjRelTol,                           &
                    AdjRCNTRL, AdjICNTRL, AdjRSTATUS, AdjISTATUS, AdjIERR)
       AdjISTATUS_T(1:6) = AdjISTATUS(1:6)

       IF (AdjIERR .LE. 0) THEN
              RETURN
       END IF

   ELSE

      !we chose to checkpoint data during the FWD integration
      !how many checkpoints were written? (at least one)
      IF (Nc .EQ. 1) THEN

          !no need to read it, everything is in the memory buffers
          !so use that information
           CALL RKINT_ADJ(N, ADY, ADJ_RHS, JACV,                    &
                    P, W, QUAD_RHS,                                 &
                    Tstart, Tend,                                   &
                    AdjAbsTol, AdjRelTol,                           &
                    AdjRCNTRL, AdjICNTRL, AdjRSTATUS, AdjISTATUS, AdjIERR)

           AdjISTATUS_T(1:6) = AdjISTATUS(1:6)

           IF (AdjIERR .LE. 0) THEN
                !an error occured, no point to continue the integration
                RETURN
           END IF

      ELSE !(Nc > 1)

          !read checkpoints in reverse order
          CALL tapes_ReadCheckpt(Nc,chk_T,chk_H,chk_Y)
          !the first time there is no need to recreate FWD trajectory for [chk_T,Tend]
          CALL RKINT_ADJ(N, ADY, ADJ_RHS, JACV,                   &
                 P, W, QUAD_RHS,                                  &
                 chk_T, Tend,                                     &
                 AdjAbsTol, AdjRelTol,                            &
                 AdjRCNTRL, AdjICNTRL, AdjRSTATUS, AdjISTATUS, AdjIERR)

          AdjISTATUS_T(1:6) = AdjISTATUS(1:6)

          IF (AdjIERR .LE. 0) THEN
                !an error occured, no point to continue the integration
                RETURN
          END IF

          chk_T_Next = chk_T

          !read the rest of the checkpoints
          !for these subintervals, the fwd trajectory needs to be recalculated
          DO idx = Nc-1,1,-1

             CALL tapes_ReadCheckpt(idx,chk_T,chk_H,chk_Y)
             CALL tapes_Flush

             !set the starting time step to allow a hot restart
             FwdRCNTRL(3) = chk_H
             FwdICNTRL(4) = 0 !not checkpointing anymore
             FwdICNTRL(5) = 1 !but still need to have memory buffering enabled

             CALL RKINT(N,chk_Y,RHS,chk_T,chk_T_Next,Nd,Nc1,          &
                    FwdAbsTol,FwdRelTol,                              &
                    FwdRCNTRL,FwdICNTRL,FwdRSTATUS,FwdISTATUS,FwdIERR)

             !accumulate statistics
             FwdISTATUS_T(1:4) = FwdISTATUS_T(1:4) + FwdISTATUS(1:4)

             IF (FwdIERR .LE. 0) THEN
                  !an error occured, no point to continue the integration
                  RETURN
             END IF

             CALL RKINT_ADJ(N, ADY, ADJ_RHS, JACV,                  &
                 P, W, QUAD_RHS,                                    &
                 chk_T, chk_T_Next,                                 &
                 AdjAbsTol, AdjRelTol,                              &
                 AdjRCNTRL, AdjICNTRL, AdjRSTATUS, AdjISTATUS, AdjIERR)

             !accumulate statistics
             AdjISTATUS_T(1:6) = AdjISTATUS_T(1:6) + AdjISTATUS(1:6)

             IF (AdjIERR .LE. 0) THEN
                  !an error occured, no point to continue the integration
                  RETURN
             END IF

             chk_T_Next = chk_T

          END DO
      END IF

   END IF

   !final update for the gradient dG/dp
   !this is needed in case the initial conditions for the forward problem depend
   !on any of the P parameters; the user must supply this routine since its
   !implementation is problem-dependent
   IF (AdjICNTRL(5) .EQ. 1) THEN
       CALL DY0DP(N,P,Tstart,ADY,Wupdate)
       CALL WAXPY(P,1.d0,Wupdate,1,W,1)
   END IF

   IF (Checkpointing) THEN
      !delete the tape files
      CALL tapes_Cleanup
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
     'Forced exit from RKINT_ADJDR due to the following error:' 

   SELECT CASE (Code)
    CASE (-1)
      PRINT * , '--> Improper checkpointing option AdjICNTRL(6)=', AdjICNTRL(6)
    CASE (-2)
      PRINT * , '--> If checkpointing is enabled, then Nc must be a positive integer.'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT

   PRINT *, "T=", T, "and H=", H

 END SUBROUTINE rkdr_ErrorMsg


 END SUBROUTINE RKINT_ADJDR


END MODULE RKADJDR_MODULE
