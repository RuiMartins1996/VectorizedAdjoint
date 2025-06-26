!Module to set the system parameters
module vanderpol_parameters
    INTEGER, parameter :: NVAR = 2          !The Van Der Pol system has two equations
    INTEGER, parameter :: NTOLS = 10        !Number of absolute/relative error tolerance values considered (it's the length of the tols array)
    INTEGER, parameter :: TREP=30           !Number of times some function is executed to obtain execution time data
    DOUBLE precision,parameter :: eps = 0.001 !The value of the \mu parameter    
    DOUBLE precision,parameter,dimension(NVAR) :: u0 = [2.0d0, -2.0d0 / 3.0d0 + 10.0d0 / (81.0 / eps) - 292.0d0 / (2187.0 / eps**2)]   !Initial conditions
    DOUBLE precision,parameter :: tols(NTOLS) = [1.0D-4, 1.0D-5,1.0D-6, 1.0D-7, 1.0D-8, 1.0D-9,1.0D-10,1.0D-11, 1.0D-12, 1.0D-13]  !Array of absolute/relative error tolerances
    DOUBLE precision,parameter :: tstart = 0.0d0,tend   = 5e-1  
end module

!Module to handle timing the algorithm
module utility
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
    subroutine  timeonce(time,n,np,tol,fun,jac,jacp,adjinit)
        use erk_adj_f90_integrator
        use vanderpol_parameters
        implicit none
        integer,intent(in) ::  n,np
        double precision,intent(out):: time

        double precision::t0,tf,tol
        external fun,jac,adjinit,jacp

        integer :: i,j
        double precision :: t1,t2,var(n),y_adj(n,n),p_adj(n,np),rtol(n),atol(n)

        do i=1,n
            var(i) = u0(i)
            rtol(i) = tol
            atol(i) = tol
        enddo

        y_adj(:,:) = 0d0;
        p_adj(:,:) = 0d0

        time = 0.0

        call cpu_time(t1)
        call integrate_adj(tin=0d0, tout=tend, NVAR=NVAR, np=np,&
            lambda=y_adj,mu = p_adj, adjinit=adjinit,y=var, rtol=rtol, atol=atol, NADJ=n,&
            fun=fun, jac=jac,jacp = jacp)
        call cpu_time(t2)

        time = t2-t1
    endsubroutine 

    !Times the execution of integrate_adj several times
    subroutine timeadjoint(n,np,tol,mean,std)
        use vanderpol_parameters
        external fun,jac,jacp,adjinit
        integer,intent(in) :: n,np
        double precision,intent(in) :: tol
        double precision,intent(out) ::mean,std
        
        integer::i
        double precision:: time,times(TREP)

        time = 0d0;

        do i=1,TREP !TREP is defined in the module vanderpol_parameters
            call timeonce(time,n,np,tol,fun,jac,jacp,adjinit)
            times(i) = time
        enddo

        call meanandstd(TREP,times,mean,std)

    endsubroutine   
end module utility





!Rhs of the Van Der Pol system
subroutine fun(n,t,x,dxdt)
    use vanderpol_parameters
    implicit none
    integer,intent(in) :: n
    double precision :: t,x(n),dxdt(n)
    
    dxdt(1) = x(2);
    dxdt(2) = ((1.0d0 - x(1)*x(1))*x(2) - x(1))/eps;
end subroutine fun

!Jacobian of fun with respect to x
subroutine jac(n, t, x, fjac)
    use vanderpol_parameters
    integer,intent(in) :: n
    double precision :: t,x(n),fjac(n,n)

    fjac(1,1) = 0.0d0;
    fjac(1,2) = 1.0d0;
    
    fjac(2,1) = - (2.0d0 * x(1) * x(2) + 1.0d0)/eps;
    fjac(2,2) = (1.0d0 - x(1) * x(1))/eps;
end subroutine jac

!Jacobian of fun with respect to \mu
subroutine jacp(n,np,t,x,fpjac)
    use vanderpol_parameters
    integer, intent(in) :: n,np
    double precision, intent(in) :: t,x(n)
    double precision, intent(inout) :: fpjac(n,np)

    fpjac(1,1) = 0.0d0;
    fpjac(2,1) = (-1.d0 / eps**2)*((1.0d0 - x(1) * x(1)) * x(2) - x(1));

end subroutine jacp

!Routine to initialize the adjoint variables lambda and mu
subroutine adjinit(n,np,NADJ,t,y,lambda,mu)
    integer :: n,np,NADJ,k
    double precision :: t,y(n),lambda(n,NADJ)
    double precision,optional :: mu(np,NADJ)
    !~~~> if number of parameters is not zero, extra adjoint varaible mu should be
    !defined
    if(NP>0 .and. .not. present(mu)) stop 'undefined argument mu'       
    !~~~>  the adjoint values at the final time
    lambda(:,:) = 0.0d0
    do k=1,NADJ
        lambda(k,k) = 1.0d0
    end do

    mu(:,:) = 0.0d0
end subroutine adjinit



!Main program
program vanderpol
    use erk_adj_f90_integrator
    use vanderpol_parameters
    use utility

    implicit none

    external fun,jac,jacp,adjinit
    integer :: i,n,k
    
    integer, parameter :: NADJ = NVAR,NPAR = 1

    double precision :: y_adj(NVAR,NADJ),var(NVAR),p_adj(NPAR,NADJ)
    double precision :: adjoints(NVAR*NPAR)
    
    double precision :: atol(NVAR),rtol(NVAR)

    double precision ::mvalues(size(tols)),stdvalues(size(tols)),mvalue,stdvalue

    integer :: ios
    character(len=128) :: fmt

    !Time the adjoint computation for all tolerances
    !------------------------------------------------
    write(*,*) "Compute execution times"
    do n=1,size(tols)
        call timeadjoint(NVAR,NPAR,tols(n),mvalue,stdvalue)
        mvalues(n) = mvalue
        stdvalues(n) = stdvalue
        write(*,*) mvalue,stdvalue
    enddo 
    write(*,*) "End computing execution times"

    !Open data file where the execution times will be written
    open(unit = 11, file = 'data.dat')  
    do i = 1,size(tols)
        write(11,*)tols(i),mvalues(i),stdvalues(i)
    enddo
    close(11)
    !------------------------------------------------


    !Loop over all tolerances and compute the solution and sensitivities for each one
    !------------------------------------------------

    !Open data file where the tolerance, solution and sensitivities will be written
    open(unit = 12, file = 'data_vanderpol_dopri5_fatode.dat')

    do n=1,size(tols)
    
        do i=1,NVAR
            rtol(i) = tols(n)
            atol(i) = tols(n)
        end do

        !Assign initial conditions
        var(:) = 0d0
        do i=1,NVAR
            var(i) = u0(i)
        end do

        p_adj(:,:) = 0d0
        y_adj(:,:) = 0d0

        call integrate_adj(tin=tstart, tout=tend, NVAR=NVAR, np=NPAR,&
            lambda=y_adj,mu = p_adj, adjinit=adjinit,y=var, rtol=rtol, atol=atol, NADJ=NADJ,&
            fun=fun, jac=jac,jacp = jacp)

        !Write the data of p_adj to the 1x2 array adjoints
        adjoints(:) = 0d0

        do i = 1,NVAR
            do k = 1,NPAR
                adjoints((i-1)*NPAR+k) = p_adj(k,i)
            end do 
        enddo
        
        ! Use ES24.16 to force scientific notation with 16 digits after the decimal
        fmt = '(ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16)'


        ! write with guaranteed full precision
        write(12, fmt, iostat=ios) tols(n), (var(i), i=1,NVAR), (adjoints(i), i=1,NVAR*NPAR)
        if (ios /= 0) then
            write(*,*) 'Error writing output, IOSTAT=', ios
        end if

        !Write the solution and sensitivities for the given tolerance
        !write(12,*) tols(n),(var(i),i=1,NVAR),(adjoints(i),i=1,NVAR*NPAR)
        
        !Print information to the console
        !write(*,*) "Solution of ODE:"
        !write(*,*) (var(i),i=1,NVAR)
        
        !write(*,*) "State Adjoints matrix"
        !do k=1,NVAR
        !    write(*,*) (y_adj(k,i),i=1,NVAR)
        !enddo 

        !write(*,*) "Parameter Adjoints matrix"
        !do k=1,NPAR
        !    write(*,*) (p_adj(k,i),i=1,NVAR)
        !enddo 
    end do

    close(12)!close data file
    !------------------------------------------------

    
end program vanderpol