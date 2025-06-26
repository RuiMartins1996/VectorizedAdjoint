############################################################################
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   Feb 27, 2011
#
#  Modified:
#
#  Send bug reports, comments or suggestions to zhang@vt.edu
#
############################################################################


export ARCH=ar -cr
export RANLIB=ranlib

# Configuration flags for linear solvers. Must choose one from the following three options.
# #  -DFULL_ALGEBRA            use BLAS and LAPACK library (full algebra)
# #  -DSPARSE_UMF              use UMFpack (sparse format) 
# #  -DSPARSE_LU               use SuperLU (sparse format)

LS_CONFIG  =  -DFULL_ALGEBRA



#~~~> gfortran (GNU FORTRAN Compiler)
export FC = gfortran
export FFLAGS = -cpp -Ofast -ffree-line-length-none $(LS_CONFIG)
BLAS = /lib/x86_64-linux-gnu/libblas.so.3
LAPACK=/lib/x86_64-linux-gnu/liblapack.so.3

FATDIR     = ../FATODE_1.2
LIBDIR     = $(FATDIR)/LSS_LIBS/

MODEL      = ADJ
FAMILY     = ERK_ADJ
APP        = vanderpol
INTEGRATOR = $(FAMILY)_f90_Integrator.o


all: $(APP)

#Compile and link
$(APP): vanderpol.o $(LSSOLVER) $(INTEGRATOR) 
	$(FC) $(FFLAGS) -o $(APP) $< $(INTEGRATOR) $(LAPACK) $(BLAS) $(XERBLA) $(LIB)

#Adjoint integrator
$(FAMILY)_f90_Integrator.o: $(FATDIR)/$(MODEL)/$(FAMILY)/$(FAMILY)_f90_Integrator.F90
	$(FC) $(FFLAGS) -c $<
	
#Compile program
vanderpol.o: vanderpol.F90 $(INTEGRATOR)
	$(FC) $(FFLAGS) -c $<

LIB = -lm

purge:  clean

clean:
	rm -f *~ *.mod *.o *.dat

help:
	@$(ECHO) "usage: make ?"














