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

MODULE UTILS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Contains implementations of BLAS-like functions for array operations.
!  Available subroutines:
!     WLAMCH   -> computes epsilon machine
!     WCOPY    -> performs vector copy operations
!     WAXPY    -> implementation of DAXPY (level 1 BLAS operation)
!     WSCAL    -> rescales a vector by a constant (similar to DSCAL)
!     SET2ZERO -> set all components of a vector to zero
!     WDOT     -> computes the scalar product of two vectors
!
!  It is advisable to replace these by an optimized BLAS implementation
!  for the target computer architecture, if such an implementation is
!  available.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  PUBLIC
  SAVE

 CONTAINS

!--------------------------------------------------------------
DOUBLE PRECISION FUNCTION WLAMCH( C )
!--------------------------------------------------------------
!  Input:
!    C (character) : any character
!  Output:
!    Epsilon machine (roundoff)
!
!  For maximum performance, it is advisable to replace this by 
!  the function from an optimized LAPACK implementation:
!          CALL SLAMCH('E') or CALL DLAMCH('E')
!--------------------------------------------------------------

      CHARACTER ::  C
      INTEGER    :: i
      DOUBLE PRECISION, SAVE  ::  Eps
      DOUBLE PRECISION  ::  Suma
      DOUBLE PRECISION, PARAMETER  ::  ONE=1.D0, HALF=0.5D0
      LOGICAL, SAVE   ::  First=.TRUE.

      IF (First) THEN
        First = .FALSE.
        Eps = HALF**(16)
        DO i = 17, 80
          Eps = Eps*HALF
          Suma = ONE + Eps
          IF (Suma.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1
      END IF

      WLAMCH = Eps

END FUNCTION WLAMCH


!--------------------------------------------------------------
SUBROUTINE WCOPY(N,X,incX,Y,incY)
!--------------------------------------------------------------
!     Copies a vector, X, to a vector, Y:  Y <- X
!     Works only for incX=incY=1
! Input:
!        N (integer) : size of the vectors X, Y
!        X (N by 1 array, double precision) : source vector
!        Y (N by 1 array, double precision) : destination vector
!        incX (integer) : storage spacing between elements of X
!                        (not used, default value of 1 used instead)
!        incY (integer) : storage spacing between elements of Y
!                        (not used, default value of 1 used instead)     
!
! Output:
!        Y (N by 1 array, double precision) : Y = X
!
!     User should replace this by the function from the optimized BLAS implementation:
!         CALL  SCOPY(N,X,1,Y,1)   or   CALL  DCOPY(N,X,1,Y,1)
!--------------------------------------------------------------

      INTEGER  :: i,incX,incY,M,MP1,N
      DOUBLE PRECISION :: X(N),Y(N)

      IF (N.LE.0) RETURN

      M = MOD(N,8)
      IF( M .NE. 0 ) THEN
        DO i = 1,M
          Y(i) = X(i)
        END DO
        IF( N .LT. 8 ) RETURN
      END IF
      MP1 = M+1
      DO i = MP1,N,8
        Y(i) = X(i)
        Y(i + 1) = X(i + 1)
        Y(i + 2) = X(i + 2)
        Y(i + 3) = X(i + 3)
        Y(i + 4) = X(i + 4)
        Y(i + 5) = X(i + 5)
        Y(i + 6) = X(i + 6)
        Y(i + 7) = X(i + 7)
      END DO

END SUBROUTINE WCOPY


!--------------------------------------------------------------
SUBROUTINE WAXPY(N,Alpha,X,incX,Y,incY)
!--------------------------------------------------------------
!     Constant times a vector plus a vector: Y <- Y + Alpha*X
!  (SAXPY-type operation)
!     Works only for incX=incY=1
!
! Input:
!        N (integer) : size of X, Y
!        X (N by 1 array, double precision) : first vector
!        Y (N by 1 array, double precision) : second vector
!        incX (integer) : storage spacing between elements of X
!                        (not used, default value of 1 used instead)
!        incY (integer) : storage spacing between elements of Y
!                        (not used, default value of 1 used instead)
!
!  Output: Y (N by 1 array, double precision) : Y = Y + Alpha*X
!
!  The user can reroute this call to the corresponding subroutine
!  in an optimized BLAS implementation, e.g.:
!
!         CALL SAXPY(N,Alpha,X,1,Y,1) or  CALL DAXPY(N,Alpha,X,1,Y,1)
!--------------------------------------------------------------

      INTEGER  :: i,incX,incY,M,MP1,N
      DOUBLE PRECISION :: X(N),Y(N),Alpha
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.d0

      IF (Alpha .EQ. ZERO) RETURN
      IF (N .LE. 0) RETURN

      M = MOD(N,4)
      IF( M .NE. 0 ) THEN
        DO i = 1,M
          Y(i) = Y(i) + Alpha*X(i)
        END DO
        IF( N .LT. 4 ) RETURN
      END IF
      MP1 = M + 1
      DO i = MP1,N,4
        Y(i) = Y(i) + Alpha*X(i)
        Y(i + 1) = Y(i + 1) + Alpha*X(i + 1)
        Y(i + 2) = Y(i + 2) + Alpha*X(i + 2)
        Y(i + 3) = Y(i + 3) + Alpha*X(i + 3)
      END DO

END SUBROUTINE WAXPY


!--------------------------------------------------------------
SUBROUTINE WSCAL(N,Alpha,X,incX)
!--------------------------------------------------------------
!     Constant times a vector (rescales the vector): x(1:N) <- Alpha*x(1:N) 
!     Works only for incX=incY=1
!
! Inputs:
!        N (integer) : size of vector X
!        Alpha (double precision) : scaling constant
!        X (N by 1 array, double precision) : vector to be rescaled
!        incX (integer) : storage spacing between elements of X
!                        (not used, default value of 1 used instead)
!
! Output:
!       X (N by 1 array, double precision) : vector rescaled by Alpha
!
!     User should eplace this with the function from an optimized BLAS implementation:
!         CALL SSCAL(N,Alpha,X,1) or  CALL DSCAL(N,Alpha,X,1)
!--------------------------------------------------------------

      INTEGER  :: i,incX,M,MP1,N
      DOUBLE PRECISION  :: X(N),Alpha
      DOUBLE PRECISION, PARAMETER  :: ZERO=0.d0, ONE=1.d0

      IF (Alpha .EQ. ONE) RETURN
      IF (N .LE. 0) RETURN

      M = MOD(N,5)
      IF( M .NE. 0 ) THEN
        IF (Alpha .EQ. (-ONE)) THEN
          DO i = 1,M
            X(i) = -X(i)
          END DO
        ELSEIF (Alpha .EQ. ZERO) THEN
          DO i = 1,M
            X(i) = ZERO
          END DO
        ELSE
          DO i = 1,M
            X(i) = Alpha*X(i)
          END DO
        END IF
        IF( N .LT. 5 ) RETURN
      END IF
      MP1 = M + 1
      IF (Alpha .EQ. (-ONE)) THEN
        DO i = MP1,N,5
          X(i)     = -X(i)
          X(i + 1) = -X(i + 1)
          X(i + 2) = -X(i + 2)
          X(i + 3) = -X(i + 3)
          X(i + 4) = -X(i + 4)
        END DO
      ELSEIF (Alpha .EQ. ZERO) THEN
        DO i = MP1,N,5
          X(i)     = ZERO
          X(i + 1) = ZERO
          X(i + 2) = ZERO
          X(i + 3) = ZERO
          X(i + 4) = ZERO
        END DO
      ELSE
        DO i = MP1,N,5
          X(i)     = Alpha*X(i)
          X(i + 1) = Alpha*X(i + 1)
          X(i + 2) = Alpha*X(i + 2)
          X(i + 3) = Alpha*X(i + 3)
          X(i + 4) = Alpha*X(i + 4)
        END DO
      END IF

END SUBROUTINE WSCAL


!--------------------------------------------------------------
SUBROUTINE SET2ZERO(N,Y)
!--------------------------------------------------------------
!    Copies zeros into the vector y:  y <- 0
!
! Input:
!       N (integer) : size of vector Y
!       Y (N by 1 array, double precision) : input vector
!
! Output:
!       Y (N by 1 array, double precision) : zero vector
!                             (i.e. Y(i)=0.d0, for all i=1..N)
!--------------------------------------------------------------

      INTEGER ::  i,M,MP1,N
      DOUBLE PRECISION, INTENT(OUT) ::  Y(N)
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0

      IF (N.LE.0) RETURN

      M = MOD(N,8)
      IF( M .NE. 0 ) THEN
        DO i = 1,M
          Y(i) = ZERO
        END DO
        IF( N .LT. 8 ) RETURN
      END IF
      MP1 = M+1
      DO i = MP1,N,8
        Y(i)     = ZERO
        Y(i + 1) = ZERO
        Y(i + 2) = ZERO
        Y(i + 3) = ZERO
        Y(i + 4) = ZERO
        Y(i + 5) = ZERO
        Y(i + 6) = ZERO
        Y(i + 7) = ZERO
      END DO

END SUBROUTINE SET2ZERO


!--------------------------------------------------------------
DOUBLE PRECISION FUNCTION WDOT(N,X,incX,Y,incY)
!--------------------------------------------------------------
!~~~>    Dot product: WDOT = x(1:N)'*y(1:N)
!~~~>    Works only for incX = incY = 1
!~~~>    Implemented after BLAS subroutines SDOT/DDOT
!
! Input:
!        N (integer) : size of X, Y
!        X (N by 1 array, double precision) : first vector
!        Y (N by 1 array, double precision) : second vector
!        incX (integer) : storage spacing between elements of X
!                        (not used, default value of 1 used instead)
!        incY (integer) : storage spacing between elements of Y
!                        (not used, default value of 1 used instead)
!
!  Output: WDOT (double precision) = x(1:N)'*y(1:N)
!
!~~~>   The user can replace this by the function from his/her
!       optimized BLAS implementation:
!          CALL SDOT(N,X,1,Y,1) or  CALL DDOT(N,X,1,Y,1)
!--------------------------------------------------------------

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, incX, incY
      DOUBLE PRECISION, INTENT(IN) :: X(N), Y(N)
      INTEGER :: M, MP1, i
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.d0

      WDOT = ZERO

      M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40

      DO i = 1,M
         WDOT = WDOT + X(i)*Y(i)
      END DO
      IF (N .LT. 5) RETURN

 40   MP1 = M + 1
      DO i = MP1,N,5
          WDOT = WDOT + X(i)*Y(i)  +  &
                    X(i+1)*Y(i+1)  +  &
                    X(i+2)*Y(i+2)  +  &
                    X(i+3)*Y(i+3)  +  &
                    X(i+4)*Y(i+4)
      END DO
      RETURN

END FUNCTION WDOT


!End of UTILS module
END MODULE UTILS