!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Module to set the system parameters
module vanderpol_parameters
  INTEGER, parameter :: NVAR = 2          !The Van Der Pol system has two equations
  INTEGER, parameter :: NPAR = 1          !The Van Der Pol system has one parameter
  INTEGER, parameter :: NTOLS = 10        !Number of absolute/relative error tolerance values considered (it's the length of the tols array)
  INTEGER, parameter :: TREP=30          !Number of times some function is executed to obtain execution time data
  DOUBLE precision,parameter :: eps = 0.001 !The value of the \mu parameter    
  DOUBLE precision,parameter,dimension(NVAR) :: u0 = [2.0d0, -2.0d0 / 3.0d0 + 10.0d0 / (81.0 / eps) - 292.0d0 / (2187.0 / eps**2)]   !Initial conditions
  DOUBLE precision,parameter :: tols(NTOLS) = [1.0D-4, 1.0D-5,1.0D-6, 1.0D-7, 1.0D-8, 1.0D-9,1.0D-10,1.0D-11, 1.0D-12, 1.0D-13]  !Array of absolute/relative error tolerances
  DOUBLE precision,parameter :: tStart = 0.0d0,tEnd   = 5e-1  
end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MODULE INTEGRATE_MODULE
IMPLICIT NONE
contains

SUBROUTINE INTEGRATE(N, P, Y, ADY, W, TIN, TOUT,TOL, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U, &
  adICNTRL_U, adRCNTRL_U, adISTATUS_U, adRSTATUS_U, adIERR_U)

   USE RK_MODULE
   USE RKADJDR_MODULE
   IMPLICIT NONE
   EXTERNAL F, ADJRHS, QUADRHS, JACEVAL_Y, DY0DP

   INTEGER :: N
   INTEGER :: P
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N), ADY(N), W(P)
   DOUBLE PRECISION, INTENT(INOUT) :: TIN  ! Start Time
   DOUBLE PRECISION, INTENT(INOUT) :: TOUT ! End Time
   DOUBLE PRECISION :: TOL

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
   Nd = 100 !checkpoint every Nd steps

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
   ICNTRL(3) = RK5
   adICNTRL(3) = RK5

!~~~> Default value: adICNTRL(4) = 0: dense output enabled
!~~~> Uncomment any of the following two lines to enable Hermite interpolation
   !adICNTRL(4) = 1 !3rd order Hermite
   !adICNTRL(4) = 2 !5th order Hermite

!~~~> Problem is parameter-dependent (eps)
   adICNTRL(5) = 1

!~~~> Enable checkpointing and memory buffering
   ICNTRL(4) = 1 !enable file checkpointing
   ICNTRL(5) = 1 !enable memory buffering
   adICNTRL(6) = 1 !inform RKINT_ADJDR that checkpointing
                   !was used in FWD mode

   ! --- REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
   RTOL(:)=tol
   ATOL(:)=RTOL(:)
   adRTOL(:)=tol
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
       TapeSize = Nd + 1
       CALL rk_AllocateTapes(CheckptEnabled,TapeSize,N)
   ELSE
       TapeSize = ICNTRL(2)
       CALL rk_AllocateTapes(CheckptEnabled,TapeSize,N)
   END IF

   !Forward integration
   CALL RKINT(N,Y,F,TIN,TOUT,  &
         Nd, Nc,                   &
         ATOL,RTOL,                &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)

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

   !IF (CheckptEnabled) THEN
   !    PRINT *, 'Wrote ', Nc, ' checkpoint(s) during FWD integration'
   !END IF
   
   !Adjoint integration
   CALL RKINT_ADJDR(N, ADY, ADJRHS, JACEVAL_Y,               &
           P, W, QUADRHS, DY0DP, F,                      &
           TIN, TOUT, Nc,                                    &
           ATOL, RTOL,                                       &
           adATOL,adRTOL,                                    &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR,               &
           adRCNTRL,adICNTRL,adRSTATUS,adISTATUS,adIERR)
   
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

! This routines computes the adjoints w.r.t the parameters for the two components of the solution
SUBROUTINE INTEGRATE_ADJ(Y,ADJOINTS, TIN, TOUT,TOL, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U, &
  adICNTRL_U, adRCNTRL_U, adISTATUS_U, adRSTATUS_U, adIERR_U)

  use vanderpol_parameters

  DOUBLE PRECISION,INTENT(OUT) :: ADJOINTS(NVAR*NPAR)
  DOUBLE PRECISION,INTENT(OUT) :: Y(NVAR)
  DOUBLE PRECISION :: ADY(NVAR), W(NPAR)
  DOUBLE PRECISION, INTENT(INOUT) :: TIN  ! Start Time
  DOUBLE PRECISION, INTENT(INOUT) :: TOUT ! End Time
  DOUBLE PRECISION :: TOL

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

  DOUBLE PRECISION :: TI,TF

  Y(1) = u0(1)
  Y(2) = u0(2)
  ADY(1) = 1.D0
  ADY(2) = 0.D0
  W = 0.d0

  ti = tStart
  tf = tEnd

  CALL INTEGRATE(NVAR, NPAR, Y, ADY, W, ti, tf,tol)  

  adjoints(1) = W(1)

  Y(1) = u0(1)
  Y(2) = u0(2)
  ADY(1) = 0.D0
  ADY(2) = 1.D0
  W = 0.d0

  ti = tStart
  tf = tEnd

  CALL INTEGRATE(NVAR, NPAR, Y, ADY, W, ti, tf,tol)

  adjoints(2) = W(1)

END SUBROUTINE INTEGRATE_ADJ

END MODULE INTEGRATE_MODULE



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module UTILITY      !Module to handle timing the algorithm
  use INTEGRATE_MODULE
  use vanderpol_parameters
  implicit none
contains
  !Computes mean value and standard deviation of the array data
  subroutine meanandstd(n,data,mean,std)
      integer :: n,i
      double precision :: mean,std
      double precision :: data(n)

      mean = 0d0
      do i = 1,n
          mean = mean+data(i)
      enddo

      mean = mean/n

      std = 0 
      do i = 1,n
          std = std+(data(i)-mean)*(data(i)-mean)
      enddo

      std = sqrt(std/(n-1))
  end subroutine meanandstd

  !Times the execution of integrate_adj once
  subroutine timeonce(time,tol)   
    use vanderpol_parameters

    DOUBLE PRECISION,INTENT(OUT) :: time
    DOUBLE PRECISION :: ti,tf  
    DOUBLE PRECISION :: TOL
    DOUBLE PRECISION :: adjoints(NPAR*NVAR)
    DOUBLE PRECISION :: y(NVAR)

    call cpu_time(ti) ! Start time
    call INTEGRATE_ADJ(y,adjoints,ti,tf,tol)
    call cpu_time(tf) ! End time

    time = tf - ti ! Calculate elapsed time
  
  endsubroutine 

  !Times the execution of integrate_adj several times
  subroutine timeadjoint(tol,mean,std)
      use vanderpol_parameters
      double precision,intent(in) :: tol
      double precision,intent(out) ::mean,std
      
      integer::i
      double precision:: time,times(TREP)

      time = 0d0;

      do i=1,TREP !TREP is defined in the module vanderpol_parameters
        call timeonce(time,tol)  
        times(i) = time
      enddo

      call meanandstd(TREP,times,mean,std)

  endsubroutine   
end module UTILITY
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Solves the Vanderpol oscillator and its first order adjoint ODE system.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PROGRAM MAIN

  use vanderpol_parameters
  use UTILITY
  use INTEGRATE_MODULE

  IMPLICIT NONE

  INTEGER :: i,j

  DOUBLE PRECISION :: ti,tf,time
  DOUBLE PRECISION :: mean,std,mvalues(NTOLS),stdvalues(NTOLS)

  DOUBLE PRECISION :: tol
  
  DOUBLE PRECISION :: adjoints(NVAR*NPAR)
  DOUBLE PRECISION :: y(NVAR)

  integer :: ios
  character(len=128) :: fmt

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Begin program
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  !Time the adjoint computation for all tolerances
  !------------------------------------------------
  write(*,*) "Compute execution times"
  do i=1,size(tols)
      tol = tols(i)
      call timeadjoint(tol,mean,std)
      mvalues(i) = mean
      stdvalues(i) = std
      write(*,*) mean,std
  enddo 
  write(*,*) "End computing execution times"

  !Write the execution times
  open(unit = 11, file = 'data.dat')  
  do i = 1,size(tols)
      write(11,*)tols(i),mvalues(i),stdvalues(i)
  enddo
  close(11)
  !------------------------------------------------

  
  !Loop over all tolerances and compute the solution and sensitivities for each one
    !------------------------------------------------

    !Open data file where the tolerance, solution and sensitivities will be written
  open(unit = 12, file = 'data_vanderpol_dopri5_denserks.dat')


  ! Use ES24.16 to force scientific notation with 16 digits after the decimal
  fmt = '(ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16)'

  write(*,*) "Computing the solution and adjoints for each tolerance"
  write(*,*) "------------------------------------------------"
  do i=1,NTOLS
  
      tol = tols(i)

      ti = tStart
      tf = tEnd
      call INTEGRATE_ADJ(y,adjoints,ti,tf,tol)
    
      !print *, 'Adjoint w.r.t. the first component of the solution:', adjoints(1)
      !print *, 'Adjoint w.r.t. the second component of the solution:', adjoints(2)

      ! write with guaranteed full precision
      write(12, fmt, iostat=ios) tols(i),(y(j),j=1,NVAR),(adjoints(j),j=1,NVAR*NPAR)
      if (ios /= 0) then
          write(*,*) 'Error writing output, IOSTAT=', ios
      end if

  end do

  close(12)!close data file
  !------------------------------------------------


END PROGRAM MAIN




!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE F(n,t,x,dxdt)
  use vanderpol_parameters
  implicit none
  integer,intent(in) :: n
  double precision :: t,x(n),dxdt(n)
  
  dxdt(1) = x(2);
  dxdt(2) = ((1.0d0 - x(1)*x(1))*x(2) - x(1))/eps;

END SUBROUTINE F


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE ADJRHS(N, T1, Y1, ADY1, ADJ1)
  use vanderpol_parameters
  IMPLICIT NONE
  integer, intent(in) :: N
  double precision, intent(in) :: t1, y1(N), ady1(N)
  double precision, intent(out) :: adj1(N)
  
  adj1(1) = -(-2.d0*y1(1)*y1(2)-1.d0)*ady1(2)/eps
  adj1(2) = -ady1(1) - (1.d0-y1(1)**2)*ady1(2)/eps
  
END SUBROUTINE ADJRHS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE JACEVAL_Y(N, T1, Y1, DY1, JAC1)
    use vanderpol_parameters
    IMPLICIT NONE
    integer, intent(in) :: N
    double precision, intent(in) :: t1, y1(N), dy1(N)
    double precision, intent(out) :: jac1(N)
  
    jac1(1) = dy1(2)
    jac1(2) = (-2.d0*y1(1)*y1(2)-1.d0)*dy1(1)/eps + (1.d0 - y1(1)**2)*dy1(2)/eps
  
END SUBROUTINE JACEVAL_Y

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE QUADRHS(N,P,T1,Y1,ADY1,W1)
  use vanderpol_parameters

  IMPLICIT NONE
  INTEGER :: N
  INTEGER :: P
  DOUBLE PRECISION, INTENT(IN) :: T1
  DOUBLE PRECISION, INTENT(IN) :: Y1(N), ADY1(N)
  DOUBLE PRECISION, INTENT(OUT) :: W1(P)
  DOUBLE PRECISION :: dfdp(2)

  !w' = -(df/dp)^T * w - (dg/dp)^T

  dfdp(1) = 0.d0
  dfdp(2) = (-1.d0 / eps**2) * ((1.d0 - y1(1)**2)*y1(2) - y1(1))

  W1 = - dfdp(1)*ADY1(1) - dfdp(2)*ADY1(2)
  
END SUBROUTINE QUADRHS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE DY0DP(N,P,T0,ADY0,GRAD0)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N, P
  DOUBLE PRECISION, INTENT(IN) :: T0
  DOUBLE PRECISION, INTENT(IN)  :: ADY0(N)
  DOUBLE PRECISION, INTENT(OUT) :: GRAD0(P)

  GRAD0 = 0.d0

END SUBROUTINE DY0DP