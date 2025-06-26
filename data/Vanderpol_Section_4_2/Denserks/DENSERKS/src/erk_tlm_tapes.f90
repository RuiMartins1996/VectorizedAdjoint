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

MODULE TLM_TAPES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Implementation of the memory-based buffers and the file checkpointing
! mechanism for the tangent linear model integration.
!
! The buffers are designed to store TLM-related information (state, stages,
! RHS values, time steps and time moments) between two consecutive file
! checkpoints. This memory and disk storage approach
! minimizes memory requirements, therefore making integrations of large-scale
! models computationally feasible.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   USE UTILS
   IMPLICIT NONE
   PUBLIC
   SAVE

!--> Units for the checkpoint files
   INTEGER, PARAMETER :: chkFileH = 205
   INTEGER, PARAMETER :: chkFileT = 206
   INTEGER, PARAMETER :: chkFileY = 207
   INTEGER, PARAMETER :: chkFileDY = 208

!--> YTape -> stores forward model trajectory
   DOUBLE PRECISION, ALLOCATABLE :: YTape(:,:)
!--> DYTape -> stores TLM trajectory
   DOUBLE PRECISION, ALLOCATABLE :: DYTape(:,:)
!--> KTape -> stores stages computed during forward integration
   DOUBLE PRECISION, ALLOCATABLE :: KTape(:,:,:)
!--> KTape -> stores stages computed during TLM integration
   DOUBLE PRECISION, ALLOCATABLE :: KtlmTape(:,:,:)
!--> TimeTape -> stores time moments corresponding to the forward
!    trajectory values in **both** YTape and DYTape
   DOUBLE PRECISION, ALLOCATABLE :: TimeTape(:)
!--> HTape -> stores the time steps taken during the forward/TLM integration
   DOUBLE PRECISION, ALLOCATABLE :: HTape(:)
!--> FTape -> stores the RHS values computed during the forward integration
   DOUBLE PRECISION, ALLOCATABLE :: FTape(:,:)
!--> FtlmTape -> stores the RHS values computed during the TLM integration
   DOUBLE PRECISION, ALLOCATABLE :: FtlmTape(:,:)

!--> stackPtr = the current stack pointer in the tape
   INTEGER :: stackPtr
   INTEGER :: stackPtrTlm
!--> backups for the stack pointer values; used for reinitialization
!    this helps in the case of multiple adjoint model integrations (see erk_soadr_M.f90)
   INTEGER :: stackPtrBkp, stackPtrTlmBkp
!--> file checkpoint "stack pointers"
   INTEGER :: chkStackPtr
!--> maximum capacity of the tapes
   INTEGER :: allocSize
!-> dimensions of the tapes and corresponding tape entries
   INTEGER :: Size1, Size2, Size3
   LOGICAL :: DoCheckpt

 CONTAINS


SUBROUTINE tlm_tapes_Allocate(CheckEnabled, Sz1, Sz2, Sz3)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Allocates all tapes.
! Must be called before any other tape-related subroutine to ensure
! tapes are correctly allocated and initialized.
! Inputs:
!       CheckEnabled - logical flag enabling/disabling checkpointing
!                      features
!       Sz1, Sz2, Sz3 : dimensions of the tapes, as follows:
!       Sz1 = the maximum number of entries in any tapes (max capacity)
!       Sz2, Sz3 = dimensions of the individual tape entries:
!             YTape = YTape(1:Sz1,1:Sz2)
!             DYTape = DYTape(1:Sz1,1:Sz2)
!             KTape = KTape(1:Sz1,1:Sz2,1:Sz3)
!             KtlmTape = KtlmTape(1:Sz1,1:Sz2,1:Sz3)
!             FTape = FTape(1:Sz1,1:Sz2)
!             FtlmTape = FtlmTape(1:Sz1,1:Sz2)
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
   stackPtrTlm = 0
   stackPtrBkp = stackPtr
   stackPtrTlmBkp = stackPtrTlm
   !chkStackPtr = 0
   allocSize = Sz1
   Size1 = Sz1
   Size2 = Sz2
   Size3 = Sz3
   DoCheckpt = CheckEnabled

   ALLOCATE (YTape(1:Sz1,1:Sz2), STAT=info)
   IF (info .NE. 0) THEN
      PRINT*,'Failed allocation of buffer Y'; STOP
   END IF
   ALLOCATE (DYTape(1:Sz1,1:Sz2), STAT=info)
   IF (info .NE. 0) THEN
      PRINT*,'Failed allocation of buffer DY'; STOP
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
   ALLOCATE (KtlmTape(1:Sz1,1:Sz2,1:Sz3), STAT=info)
   IF (info .NE. 0) THEN
      PRINT*,'Failed allocation of buffer Ktlm'; STOP
   END IF
   ALLOCATE (FTape(1:Sz1,1:Sz2), STAT=info)
   IF (info .NE. 0) THEN
      PRINT*,'Failed allocation of buffer F'; STOP
   END IF
   ALLOCATE (FtlmTape(1:Sz1,1:Sz2), STAT=info)
   IF (info .NE. 0) THEN
      PRINT*,'Failed allocation of buffer Ftlm'; STOP
   END IF

   YTape = 0.d0
   DYTape = 0.d0
   FTape = 0.d0
   FtlmTape = 0.d0
   HTape = 0.d0
   KTape = 0.d0
   KtlmTape = 0.d0
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
     inquire(iolength=length) DYTape(1,:)
     open(chkFileDY, err=99, file='chk_DY.dat', &
             & status='replace', access='direct', &
             & form='unformatted', recl=length)
     inquire(iolength=length) TimeTape(1)
     open(chkFileT, err=99, file='chk_T.dat', &
             & status='replace', access='direct', &
             & form='unformatted', recl=length)

     chkStackPtr = 0

   END IF

   RETURN

99 PRINT *,'Could not open the tape files for writing.', & 
          & 'Make sure you have disk write access.'
   STOP 0


END SUBROUTINE tlm_tapes_Allocate



SUBROUTINE tlm_tapes_Deallocate
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

   IF (allocated(DYTape)) THEN 
       DEALLOCATE(DYTape, STAT=info)
       IF (info .NE. 0) THEN
              PRINT*,'Failed deallocation of buffer DY'; STOP
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

   IF (allocated(KtlmTape)) THEN
      DEALLOCATE(KtlmTape, STAT=info)
      IF (info .NE. 0) THEN
         PRINT*,'Failed deallocation of buffer Ktlm'; STOP
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

   IF (allocated(FtlmTape)) THEN
      DEALLOCATE(FtlmTape, STAT=info)
      IF (info .NE. 0) THEN
         PRINT*,'Failed deallocation of buffer Ftlm'; STOP
      END IF
   END IF

   stackPtr = 0
   stackPtrTlm = 0
   stackPtrBkp = 0
   stackPtrTlmBkp = 0
   allocSize = 0
   Size1 = 0
   Size2 = 0
   Size3 = 0

   IF (DoCheckpt) THEN
      !close the files
      close(chkFileH)
      close(chkFileY)
      close(chkFileDY)
      close(chkFileT)
   END IF

END SUBROUTINE tlm_tapes_Deallocate



SUBROUTINE tlm_tapes_Cleanup
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Performs final cleanup: removes the checkpoint files
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

! Can uncomment this on Linux/Unix systems with certain compilers
! Otherwise one will have to perform the cleanup of the tape
! files by hand.
!  CALL system('rm chk_Y.dat')
!  CALL system('rm chk_DY.dat')
!  CALL system('rm chk_T.dat')
!  CALL system('rm chk_H.dat')

END SUBROUTINE tlm_tapes_Cleanup




SUBROUTINE tlm_tapes_CheckAllocated(IERR)
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
   INTEGER :: reclen
   LOGICAL :: isOpen, exists

   !check if the tapes have been allocated 
   !(user has to take care of this in the driver by calling dopri_AllocateTapes)
   IF ((.NOT. allocated(YTape)) .OR. &
        & (.NOT. allocated(DYtape)) .OR. &
        & (.NOT. allocated(KTape)) .OR. &
        & (.NOT. allocated(KtlmTape)) .OR. &
        & (.NOT. allocated(HTape)) .OR. & 
        & (.NOT. allocated(TimeTape)) .OR. &
        & (.NOT. allocated(FtlmTape)) .OR. &
        & (.NOT. allocated(FTape))) THEN
           CALL tlm_tapes_ErrorMsg(-1,IERR)
           RETURN
   END IF

   IF (DoCheckpt) THEN

     !check that each file exists & if it's not open, open it
     inquire(unit=chkFileH,exist=exists,opened=isOpen)
     IF (.NOT. exists) THEN
         CALL tlm_tapes_ErrorMsg(-2,IERR)
         RETURN
     ELSEIF (.NOT. isOpen) THEN
         inquire(iolength=reclen) HTape(1)
         open(unit=chkFileH, status='old', err=98, access='direct', recl=reclen)
     END IF

     inquire(unit=chkFileY,exist=exists,opened=isOpen)
     IF (.NOT. exists) THEN
         CALL tlm_tapes_ErrorMsg(-2,IERR)
         RETURN
     ELSEIF (.NOT. isOpen) THEN
         inquire(iolength=reclen) YTape(1,:)
         open(unit=chkFileY, status='old', err=98, access='direct', recl=reclen)
     END IF

     inquire(unit=chkFileT,exist=exists,opened=isOpen)
     IF (.NOT. exists) THEN
         CALL tlm_tapes_ErrorMsg(-2,IERR)
         RETURN
     ELSEIF (.NOT. isOpen) THEN
         inquire(iolength=reclen) TimeTape(1)
         open(unit=chkFileT, status='old', err=98, access='direct', recl=reclen)
     END IF

     inquire(unit=chkFileDY,exist=exists,opened=isOpen)
     IF (.NOT. exists) THEN
         CALL tlm_tapes_ErrorMsg(-2,IERR)
         RETURN
     ELSEIF (.NOT. isOpen) THEN
         inquire(iolength=reclen) DYTape(1,:)
         open(unit=chkFileDY, status='old', err=98, access='direct', recl=reclen)
     END IF

   END IF

   IERR = 1 !successful verification

   RETURN

98 CALL tlm_tapes_ErrorMsg(-3,IERR)


END SUBROUTINE tlm_tapes_CheckAllocated



SUBROUTINE tlm_tapes_Flush
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Resets the stack pointer for the memory buffers and 
! zeroes out all buffer elements
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

  !flushes ONLY the memory-based tapes
  !leaves the checkpoint files and chkStackPtr unchanged
  stackPtr = 0
  StackPtrTlm = 0

  YTape = 0.d0
  DYTape = 0.d0
  KTape = 0.d0
  KtlmTape = 0.d0
  HTape = 0.d0
  FTape = 0.d0
  FTlmTape = 0.d0
  TimeTape = 0.d0

END SUBROUTINE tlm_tapes_Flush



SUBROUTINE tlm_tapes_ReinitPointers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Reinitializes the stack pointers for the memory buffers
! This is useful when multiple adjoint model integrations
! are performed that reuse the same forward model and
! TLM trajectory information.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

   stackPtr = stackPtrBkp
   stackPtrTlm = stackPtrTlmBkp


END SUBROUTINE tlm_tapes_ReinitPointers


SUBROUTINE tlm_tapes_WriteCheckpt(T,H,Y,DY)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Writes checkpoint data corresponding to a single checkpoint to
! the appropriate files. 
!
! Enough data is stored to allow a hot restart of the integration.
! Every entry can be viewed as a quadruple chk = (T,H,Y,DY)
!
! Inputs:
!      T (double precision) : time moment
!      H (double precision) : time step
!      Y (Size2 by 1 array, double precision) : Forward model numerical solution at T
!      DY (Size2 by 1 array, double precision) : TLM numerical solution at T
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: T
  DOUBLE PRECISION, INTENT(IN) :: H
  DOUBLE PRECISION, INTENT(IN) :: Y(Size2), DY(Size2)

  chkStackPtr = chkStackPtr + 1
  !write entry to checkpoint files
  write(chkFileH, rec=chkStackPtr) H
  write(chkFileT, rec=chkStackPtr) T
  write(chkFileY, rec=chkStackPtr) Y
  write(chkFileDY, rec=chkStackPtr) DY

END SUBROUTINE tlm_tapes_WriteCheckpt



SUBROUTINE tlm_tapes_ReadCheckpt(ChkIndx,T,H,Y,DY)
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
!      DY (Size2 by 1 array, double precision) : TLM numerical solution at T
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ChkIndx
  DOUBLE PRECISION, INTENT(OUT) :: T, H
  DOUBLE PRECISION, INTENT(OUT) :: Y(Size2), DY(Size2)

  read(chkFileH, rec=ChkIndx) H
  read(chkFileY, rec=ChkIndx) Y
  read(chkFileDY, rec=ChkIndx) DY
  read(chkFileT, rec=ChkIndx) T

END SUBROUTINE tlm_tapes_ReadCheckpt



SUBROUTINE tlm_tapes_Write(T,H,K,Ktlm,F,Ftlm,Y,DY)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Writes all the arguments in the corresponding memory buffers
! Inputs: 
!      T (double precision) : time moment
!      H (double precision) : time step
!      K (Size2 by Size3 array, double precision) : Runge-Kutta stages for forward model
!      Ktlm (Size2 by Size3 array, double precision) : Runge-Kutta stages for TLM
!      F (Size2 by 1 array, double precision) : Forward model ODE RHS at (T,Y)
!      Ftlm (Size2 by 1 array, double precision) : TLM ODE RHS at (T,DY)
!      Y (Size2 by 1 array, double precision) : Forward model numerical solution at T
!      DY (Size2 by 1 array, double precision) : TLM numerical solution at T
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: H, T
  DOUBLE PRECISION, INTENT(IN) :: K(Size2,Size3), Y(Size2), F(Size2)
  DOUBLE PRECISION, INTENT(IN) :: Ktlm(Size2,Size3), DY(Size2), Ftlm(Size2)
  INTEGER i

  !Guard against buffer overflow
  IF (stackPtr .GE. allocSize) THEN
    PRINT *,'Cannot write entry in tape. Buffer overflow! stackPtr = ', stackPtr
    STOP 1
  END IF

  stackPtr = stackPtr + 1
  stackPtrTlm = stackPtrTlm + 1
  stackPtrBkp = stackPtr
  stackPtrTlmBkp = stackPtrTlm

  HTape(stackPtr) = H
  TimeTape(stackPtr) = T

  DO i=1,Size3
     CALL WCOPY(Size2,K(:,i),1,KTape(stackPtr,:,i),1)
     CALL WCOPY(Size2,Ktlm(:,i),1,KtlmTape(stackPtr,:,i),1)
  END DO

  CALL WCOPY(Size2,Y,1,YTape(stackPtr,:),1)
  CALL WCOPY(Size2,DY,1,DYTape(stackPtr,:),1)
  CALL WCOPY(Size2,F,1,FTape(stackPtr,:),1)
  CALL WCOPY(Size2,Ftlm,1,FtlmTape(stackPtr,:),1)

END SUBROUTINE tlm_tapes_Write



SUBROUTINE tlm_tapes_LookupDense(Tau,Tc,Hc,Kc,Yc,Ycp1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Performs a lookup in the forward model buffers and returns the entries for largest Tc <= Tau
! This data is useful for implementing dense output for interpolation of the forward model solution.
!
! Input:
!       Tau (double precision) : interpolation time
! Output:
!       Tc (double precision) : largest time in the time memory buffer s.t. Tc <= Tau
!       Hc (double precision) : time step taken at Tc
!       Kc (Size2 by Size3 array, double precision) : Runge-Kutta stages (forward model)
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
          IF (PRESENT(Ycp1)) THEN
              CALL WCOPY(Size2,YTape(i+1,:),1,Ycp1,1)
          END IF

          !discard the values with TimeTape(i) > Tau
          stackPtr = i+1
          EXIT

     END IF

  END DO

END SUBROUTINE tlm_tapes_LookupDense



SUBROUTINE tlm_tapes_LookupDenseTLM(Tau,Tc,Hc,Kc,Yc,Ycp1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Performs a lookup in the TLM buffers and returns the entries for largest Tc <= Tau
! This data is useful for implementing dense output for interpolation of the forward model solution.
!
! Input:
!       Tau (double precision) : interpolation time
! Output:
!       Tc (double precision) : largest time in the time memory buffer s.t. Tc <= Tau
!       Hc (double precision) : time step taken at Tc
!       Kc (Size2 by Size3 array, double precision) : Runge-Kutta stages (tangent linear model)
!       Yc (Size2 by 1 array, double precision) : tangent linear model solution at Tc
!       Ycp1 (Size2 by 1 array, double precision, optional) : tangent linear model solution at Tc+Hc
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: Tau
  DOUBLE PRECISION, INTENT(OUT) :: Tc, Hc
  DOUBLE PRECISION, INTENT(OUT) :: Kc(Size2,Size3), Yc(Size2)
  DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: Ycp1(Size2)
  INTEGER :: i,j

  DO i= stackPtrTlm-1,1,-1

     IF (Tau >= TimeTape(i)) THEN

          Tc = TimeTape(i)
          Hc = HTape(i)
          CALL WCOPY(Size2,DYTape(i,:),1,Yc,1)
          DO j=1,Size3
                CALL WCOPY(Size2,KtlmTape(i,:,j),1,Kc(:,j),1)
          END DO
          IF (PRESENT(Ycp1)) THEN
              CALL WCOPY(Size2,DYTape(i+1,:),1,Ycp1,1)
          END IF

          !discard the values with TimeTape(i) > Tau
          stackPtrTlm = i+1
          EXIT

     END IF

  END DO

END SUBROUTINE tlm_tapes_LookupDenseTLM



SUBROUTINE tlm_tapes_LookupHermite(Tau,Tc,Tcp1,Yc,Ycp1,Fc,Fcp1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Performs a buffer lookup.
! Returns data required for 3rd or 5th-order Hermite interpolation of the forward solution.
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

END SUBROUTINE tlm_tapes_LookupHermite



SUBROUTINE tlm_tapes_LookupHermiteTLM(Tau,Tc,Tcp1,Yc,Ycp1,Fc,Fcp1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Performs a TLM buffer lookup.
! Returns data required for 3rd or 5th-order Hermite interpolation of the TLM solution.
! Input:
!       Tau (double precision) : interpolation time
! Output:
!       Tc (double precision) : largest time in the time memory buffer s.t. Tc <= Tau
!       Tcp1 (double precision) : next time moment after Tc
!       Yc (Size2 by 1 array, double precision) : TLM solution at Tc
!       Ycp1 (Size2 by 1 array, double precision, optional) : TLM solution at Tcp1
!       Fc (Size2 by 1 array, double precision) : TLM ODE RHS at Tc
!       Fcp1 (Size2 by 1 array, double precision, optional) : TLM ODE RHS at Tcp1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: Tau
  DOUBLE PRECISION, INTENT(OUT) :: Tc, Tcp1
  DOUBLE PRECISION, INTENT(OUT) :: Yc(Size2), Ycp1(Size2)
  DOUBLE PRECISION, INTENT(OUT) :: Fc(Size2), Fcp1(Size2)
  INTEGER i

  DO i= stackPtrTlm-1,1,-1

     IF (Tau >= TimeTape(i)) THEN

          Tc = TimeTape(i)
          Tcp1 = TimeTape(i+1)

          CALL WCOPY(Size2,DYTape(i,:),1,Yc,1)
          CALL WCOPY(Size2,DYTape(i+1,:),1,Ycp1,1)
          CALL WCOPY(Size2,FtlmTape(i,:),1,Fc,1)
          CALL WCOPY(Size2,FtlmTape(i+1,:),1,Fcp1,1)

          !discard the values with TimeTape(i) > Tau
          stackPtrTlm = i+1
          EXIT

     END IF

  END DO

END SUBROUTINE tlm_tapes_LookupHermiteTLM


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE tlm_tapes_ErrorMsg(Code,IERR)
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

END SUBROUTINE tlm_tapes_ErrorMsg


!End of TLM_TAPES module
END MODULE TLM_TAPES