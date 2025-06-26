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

MODULE TAPES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Implementation of the memory-based buffers (a.k.a. tapes) and the file checkpointing
! mechanism for the forward model integration.
!
! The memory buffers are designed to store forward model-related information
! (state, stages, forward ODE RHS values, time steps and time moments) between
! two consecutive file checkpoints. This memory and disk storage approach
! minimizes memory requirements, therefore making integrations of large-scale
! models computationally feasible.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   USE UTILS
   IMPLICIT NONE
   PUBLIC
   SAVE

   !units for the checkpoint files
   INTEGER, PARAMETER :: chkFileH = 201
   INTEGER, PARAMETER :: chkFileT = 202
   INTEGER, PARAMETER :: chkFileY = 203

!--> YTape -> stores forward trajectory
   DOUBLE PRECISION, ALLOCATABLE :: YTape(:,:)
!--> KTape -> stores stages computed during forward integration
   DOUBLE PRECISION, ALLOCATABLE :: KTape(:,:,:)
!--> TimeTape -> stores time moments corresponding to the forward 
   !trajectory values in YTape
   DOUBLE PRECISION, ALLOCATABLE :: TimeTape(:)
!--> HTape -> stores the time steps taken during the forward integration
   DOUBLE PRECISION, ALLOCATABLE :: HTape(:)
!--> FTape -> stores the RHS values computed during the forward integration
   DOUBLE PRECISION, ALLOCATABLE :: FTape(:,:)

!--> stackPtr = the current stack pointer in the tape
!--> chkStackPtr = stack pointer for the checkpoints; useful when writing
!    and/or reading a checkpoint based on its index
   INTEGER :: stackPtr, chkStackPtr
   INTEGER :: stackPtrBkp
!--> maximum capacity of the tapes
   INTEGER :: allocSize
!--> dimensions of the tapes and corresponding tape entries
   INTEGER :: Size1, Size2, Size3
!--> is checkpointing enabled?
   LOGICAL :: DoCheckpt

 CONTAINS

SUBROUTINE tapes_Allocate(CheckEnabled,Sz1, Sz2, Sz3)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Allocates all tapes.
! Must be called before any other tape-related subroutine to ensure
! tapes are correctly allocated and initialized.
! Input:
!       CheckEnabled = flag signaling that checkpointing is enabled
!                      and checkpoint files should be opened for I/O
!       Sz1, Sz2, Sz3 : dimensions of the tapes, as follows:
!       Sz1 = the maximum number of entries in any tapes (max capacity)
!       Sz2, Sz3 = dimensions of the individual tape entries:
!             YTape = YTape(1:Sz1,1:Sz2)
!             KTape = Ktape(1:Sz1,1:Sz2,1:Sz3)
!             FTape = FTape(1:Sz1,1:Sz2)
!             HTape = HTape(1:Sz1)
!             TimeTape = TimeTape(1:Sz1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: Sz1, Sz2, Sz3
   LOGICAL, INTENT(IN) :: CheckEnabled
   INTEGER info, length

   IF ((Sz1 .LE. 0) & 
          &  .OR. (Sz2 .LE. 0) &
          &  .OR. (Sz3 .le. 0)) THEN
      PRINT *,'User specified invalid sizes for the tape array!'
      STOP 1
   END IF

   stackPtr = 0
   chkStackPtr = 0
   allocSize = Sz1
   Size1 = Sz1
   Size2 = Sz2
   Size3 = Sz3
   DoCheckpt = CheckEnabled

   ALLOCATE (YTape(1:Sz1,1:Sz2), STAT=info)
   IF (info .NE. 0) THEN
      PRINT*,'Failed allocation of buffer Y'; STOP
   END IF
   ALLOCATE (HTape(1:Sz1), STAT=info)
   IF (info .NE. 0) THEN
      PRINT*,'Failed allocation of buffer H'; STOP
   END IF
   ALLOCATE (TimeTape(1:Sz1), STAT=info)
   IF (info .NE. 0) THEN
      PRINT*,'Failed allocation of buffer T'; STOP
   END IF
   ALLOCATE (KTape(1:Sz1,1:Sz2,1:Sz3), STAT=info)
   IF (info .NE. 0) THEN
      PRINT*,'Failed allocation of buffer K'; STOP
   END IF
   ALLOCATE (FTape(1:Sz1,1:Sz2), STAT=info)
   IF (info .NE. 0) THEN
      PRINT*,'Failed allocation of buffer F'; STOP
   END IF

   YTape = 0.d0
   FTape = 0.d0
   HTape = 0.d0
   KTape = 0.d0
   TimeTape = 0.d0

   IF (DoCheckpt) THEN
     !prepare the tape files for writes/reads
     inquire(iolength=length) HTape(1)
     open(chkFileH, err=99, file='chk_H.dat', &
             & status='replace', access='direct', &
             & form='unformatted', recl=length)
     inquire(iolength=length) YTape(1,:)
     open(chkFileY, err=99, file='chk_Y.dat', &
             & status='replace', access='direct', &
             & form='unformatted', recl=length)
     inquire(iolength=length) TimeTape(1)
     open(chkFileT, err=99, file='chk_T.dat', &
             & status='replace', access='direct', &
             & form='unformatted', recl=length)
   END IF

   RETURN

99 PRINT *,'Could not open the tape files for writing.', & 
          & 'Make sure you have disk write access.'
   STOP 0

END SUBROUTINE tapes_Allocate



SUBROUTINE tapes_Deallocate
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Deallocates all memory tapes and closes the checkpoint files.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   INTEGER info

   IF (allocated(YTape)) THEN 
       DEALLOCATE(YTape, STAT=info)
       IF (info .NE. 0) THEN
         PRINT*,'Failed deallocation of buffer Y'; STOP
       END IF
   END IF

   IF (allocated(HTape)) THEN
      DEALLOCATE(HTape, STAT=info)
      IF (info .NE. 0) THEN
         PRINT*,'Failed deallocation of buffer H'; STOP
      END IF
   END IF

   IF (allocated(KTape)) THEN
      DEALLOCATE(KTape, STAT=info)
      IF (info .NE. 0) THEN
         PRINT*,'Failed deallocation of buffer K'; STOP
      END IF
   END IF

   IF (allocated(TimeTape)) THEN
      DEALLOCATE(TimeTape, STAT=info)
      IF (info .NE. 0) THEN
         PRINT*,'Failed deallocation of buffer T'; STOP
      END IF
   END IF

   IF (allocated(FTape)) THEN
      DEALLOCATE(FTape, STAT=info)
      IF (info .NE. 0) THEN
         PRINT*,'Failed deallocation of buffer F'; STOP
      END IF
   END IF

   IF (DoCheckpt) THEN
      !close the files
      close(chkFileH)
      close(chkFileY)
      close(chkFileT)
   END IF

   stackPtr = 0
   stackPtrBkp = 0
   allocSize = 0
   chkStackPtr = 0
   Size1 = 0
   Size2 = 0
   Size3 = 0

END SUBROUTINE tapes_Deallocate



SUBROUTINE tapes_Cleanup
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Performs final cleanup: removes the checkpoint files
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  IMPLICIT NONE

! Can uncomment this on Linux/Unix systems with certain compilers
! Otherwise one will have to perform the cleanup of the tape
! files by hand.
!  CALL system('rm chk_Y.dat')
!  CALL system('rm chk_T.dat')
!  CALL system('rm chk_H.dat')

END SUBROUTINE tapes_Cleanup



SUBROUTINE tapes_ReinitPointers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Reinitializes the stack pointer for the memory buffers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

   stackPtr = stackPtrBkp

END SUBROUTINE tapes_ReinitPointers 


SUBROUTINE tapes_CheckAllocated(IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Checks if tapes have been allocated, if not signals an error
! and sets IERR (error indicator) to a nonzero value.
! It is advisable to call this function before performing any tape
! operations, to avoid any memory-related errors (e.g. buffer overflows)
!
! Output: 
!      IERR (integer) : Error code
!                 1 : no error occured (tapes and files have been properly
!                     initialized)
!                 < 0 : an error occured; IERR = error code.
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   INTEGER, INTENT(INOUT) :: IERR
   LOGICAL :: exists
   LOGICAL :: isOpen
   INTEGER :: reclen

   exists = .FALSE.
   isOpen = .FALSE.

   !check if the tapes have been allocated 
   !(user has to take care of this in the driver by calling rk_AllocateTapes)
   IF ((.NOT. allocated(YTape)) .OR. &
        & (.NOT. allocated(Ktape)) .OR. &
        & (.NOT. allocated(HTape)) .OR. & 
        & (.NOT. allocated(TimeTape)) .OR. &
        & (.NOT. allocated(FTape))) THEN
	   !call error handler subroutine
           CALL tapes_ErrorMsg(-1,IERR)
           RETURN
   END IF

   IF (DoCheckpt) THEN

      !check that each file exists; if it is not open, open it
      inquire(unit=chkFileH,exist=exists,opened=isOpen)
      IF (.NOT. exists) THEN
         CALL tapes_ErrorMsg(-2,IERR)
         RETURN
	 ELSEIF (.NOT. isOpen) THEN
         inquire(iolength=reclen) HTape(1)
         open(unit=chkFileH, status='old', err=98, access='direct', recl=reclen)
      END IF

      inquire(unit=chkFileY,exist=exists,opened=isOpen)
      IF (.NOT. exists) THEN
         CALL tapes_ErrorMsg(-2,IERR)
         RETURN
      ELSEIF (.NOT. isOpen) THEN
         inquire(iolength=reclen) YTape(1,:)
         open(unit=chkFileY, status='old', err=98, access='direct', recl=reclen)
      END IF

      inquire(unit=chkFileT,exist=exists,opened=isOpen)
      IF (.NOT. exists) THEN
         CALL tapes_ErrorMsg(-2,IERR)
         RETURN
      ELSEIF (.NOT. isOpen) THEN
         inquire(iolength=reclen) TimeTape(1)
         open(unit=chkFileT, status='old', err=98, access='direct', recl=reclen)
      END IF

   END IF

   IERR = 1 !successful verification

   RETURN

98 CALL tapes_ErrorMsg(-3,IERR)

END SUBROUTINE tapes_CheckAllocated



SUBROUTINE tapes_Flush
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Resets the stack pointer for the memory buffers and zeroes out all buffer elements
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !flushes ONLY the memory buffers
  !leaves the checkpoint files and chkStackPtr unchanged
  stackPtr = 0
  stackPtrBkp = 0

  YTape = 0.d0
  KTape = 0.d0
  HTape = 0.d0
  FTape = 0.d0
  TimeTape = 0.d0

END SUBROUTINE tapes_Flush



SUBROUTINE tapes_Write(T,H,K,F,Y)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Writes all the arguments in the corresponding memory buffers
! Inputs: 
!      T (double precision) : time moment
!      H (double precision) : time step
!      K (Size2 by Size3 array, double precision) : Runge-Kutta stages
!      F (Size2 by 1 array, double precision) : Forward model ODE RHS
!      Y (Size2 by 1 array, double precision) : Forward model numerical solution at T
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: H, T
  DOUBLE PRECISION, INTENT(IN) :: K(Size2,Size3), Y(Size2), F(Size2)
  INTEGER i

  !guard against buffer overflow
  IF (stackPtr .GE. allocSize) THEN
    PRINT *,'Cannot write entry in the memory buffers. Buffer overflow! stackPtr = ', stackPtr
    STOP 1
  END IF

  stackPtr = stackPtr + 1
  stackPtrBkp = stackPtr

  HTape(stackPtr) = H
  TimeTape(stackPtr) = T

  DO i=1,Size3
     CALL WCOPY(Size2,K(:,i),1,KTape(stackPtr,:,i),1)
  END DO

  CALL WCOPY(Size2,Y,1,YTape(stackPtr,:),1)
  CALL WCOPY(Size2,F,1,FTape(stackPtr,:),1)

END SUBROUTINE tapes_Write



SUBROUTINE tapes_WriteCheckpt(T,H,Y)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Writes checkpoint data corresponding to a single checkpoint to
! the appropriate files.
!
! Enough data is stored to allow a hot restart of the integration.
! Every entry can be viewed as a triplet chk = (T,H,Y)
!
! Inputs:
!      T (double precision) : time moment
!      H (double precision) : time step
!      Y (Size2 by 1 array, double precision) : Forward model numerical solution at T
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: T
  DOUBLE PRECISION, INTENT(IN) :: H
  DOUBLE PRECISION, INTENT(IN) :: Y(Size2)

  chkStackPtr = chkStackPtr + 1
  !write entry to checkpoint files
  write(chkFileH, rec=chkStackPtr) H
  write(chkFileT, rec=chkStackPtr) T
  write(chkFileY, rec=chkStackPtr) Y

END SUBROUTINE tapes_WriteCheckpt



SUBROUTINE tapes_ReadCheckpt(ChkIndx,T,H,Y)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Reads checkpoint data corresponding to the ChkIndx-th checkpoint.
!
! Inputs:
!      ChkIndx (integer) : index of the checkpoint that must be read
!                          (we employ 1-based checkpoint indexing)
!
! Outputs:
!      T (double precision) : time
!      H (double precision) : time step
!      Y (Size2 by 1 array, double precision) : Forward model numerical solution at T
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ChkIndx
  DOUBLE PRECISION, INTENT(OUT) :: T, H
  DOUBLE PRECISION, INTENT(OUT) :: Y(Size2)

  read(chkFileH, rec=ChkIndx) H
  read(chkFileY, rec=ChkIndx) Y
  read(chkFileT, rec=ChkIndx) T

END SUBROUTINE tapes_ReadCheckpt



SUBROUTINE tapes_LookupDense(Tau,Tc,Hc,Kc,Yc,Ycp1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Performs a buffer lookup and returns the tape entries for largest Tc <= Tau
! This data is useful for implementing dense output.
!
! Input:
!       Tau (double precision) : interpolation time
! Output:
!       Tc (double precision) : largest time in the time memory buffer s.t. Tc <= Tau
!       Hc (double precision) : time step taken at Tc
!       Kc (Size2 by Size3 array, double precision) : Runge-Kutta stages
!       Yc (Size2 by 1 array, double precision) : forward model solution at Tc
!       Ycp1 (Size2 by 1 array, double precision, optional) : Forward model solution at Tc+Hc
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: Tau
  DOUBLE PRECISION, INTENT(OUT) :: Tc, Hc
  DOUBLE PRECISION, INTENT(OUT) :: Kc(Size2,Size3), Yc(Size2)
  DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: Ycp1(Size2)
  INTEGER :: i,j

  DO i= stackPtr-1,1,-1
     IF (Tau >= TimeTape(i)) THEN
          Tc = TimeTape(i)
          Hc = HTape(i)
          CALL WCOPY(Size2,YTape(i,:),1,Yc,1)
          DO j=1,Size3
                CALL WCOPY(Size2,KTape(i,:,j),1,Kc(:,j),1)
          END DO

          !Ycp1 is needed for 7th order dense output (with DOPRI8(6))
          IF (PRESENT(Ycp1)) THEN
              CALL WCOPY(Size2,YTape(i+1,:),1,Ycp1,1)
          END IF

          !discard the values with TimeTape(i) > Tau
          stackPtr = i+1
          EXIT
     END IF
  END DO

END SUBROUTINE tapes_LookupDense



SUBROUTINE tapes_LookupHermite(Tau,Tc,Tcp1,Yc,Ycp1,Fc,Fcp1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Performs a buffer lookup.
! Returns data required for 3rd or 5th-order Hermite interpolation.
! Input:
!       Tau (double precision) : interpolation time
! Output:
!       Tc (double precision) : largest time in the time memory buffer s.t. Tc <= Tau
!       Tcp1 (double precision) : next time moment after Tc
!       Yc (Size2 by 1 array, double precision) : forward model solution at Tc
!       Ycp1 (Size2 by 1 array, double precision, optional) : Forward model solution at Tcp1
!       Fc (Size2 by 1 array, double precision) : forward model ODE RHS at Tc
!       Fcp1 (Size2 by 1 array, double precision, optional) : Forward model ODE RHS at Tcp1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: Tau
  DOUBLE PRECISION, INTENT(OUT) :: Tc, Tcp1
  DOUBLE PRECISION, INTENT(OUT) :: Yc(Size2), Ycp1(Size2)
  DOUBLE PRECISION, INTENT(OUT) :: Fc(Size2), Fcp1(Size2)
  INTEGER i

  DO i= stackPtr-1,1,-1
     IF (Tau >= TimeTape(i)) THEN

          Tc = TimeTape(i)
          Tcp1 = TimeTape(i+1)

          CALL WCOPY(Size2,YTape(i,:),1,Yc,1)
          CALL WCOPY(Size2,YTape(i+1,:),1,Ycp1,1)
          CALL WCOPY(Size2,FTape(i,:),1,Fc,1)
          CALL WCOPY(Size2,FTape(i+1,:),1,Fcp1,1)

          !discard the values with TimeTape(i) > Tau
          stackPtr = i+1
          EXIT
     END IF
  END DO

END SUBROUTINE tapes_LookupHermite


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE tapes_ErrorMsg(Code,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all tape-related error messages
!
! Input:
!      Code (integer) : error code
! Output:
!      IERR (integer) : Error indicator. Upon exit from tapes_ErrorMsg, IERR = Code.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR

   IERR = Code

   SELECT CASE (Code)

      CASE (-1)
        PRINT * , '--> Tapes have not been allocated.'
      CASE (-2)
        PRINT * , '--> One or more checkpoint files do not exist.'
      CASE (-3)
        PRINT * , '--> Could not open one or more checkpoint files. ', &
                  & 'Make sure you have marked the files for read access.'
      CASE DEFAULT
        PRINT *, 'Unknown Error code: ', Code

   END SELECT

END SUBROUTINE tapes_ErrorMsg


!End of TAPES module
END MODULE TAPES