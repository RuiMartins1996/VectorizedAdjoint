############################################################################
#
#  Program:         DENSERKS
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   August 28, 2007
#
#  Modified:
#
#  Send bug reports, comments or suggestions to mihai@vt.edu
#
############################################################################

# Use ar unless compiling with ifort.
# If using ifort, enable xiar (see below)
export ARCH=ar -cr
export RANLIB=ranlib

#~~~> Intel Fortran compiler ()
#~~~> We maintain floating point precision using -mp
#     This disables some arithmetic optimizations, but ensures corectness of the resulting binary code
#export FC=ifort
#export FFLAGS=-fast -mp -msse3 -vec-report0
#export MODSWITCH=-module
#~~~> Use xiar instead of ar when compiling with -ipo (-fast implicitly enables -ipo)
#     This is for use with the Intel Fortran compiler (ifort) <only>. 
#     With other compilers, one should use the usual ar archival tool.
#export ARCH=xiar -qipo -cr

#~~~> Linker options for Intel Fortran
#     Try this in case ifort has problem with the built-in mathematical functions
# export LOPT=-lsvml


export FC=gfortran
#export FFLAGS= -Ofast -ffree-line-length-none $(LS_CONFIG)
export FFLAGS= -g -ffree-line-length-none $(LS_CONFIG)

export MODSWITCH=-I

#path to the DENSERKS library
DENSERKLIB=../DENSERKS/lib/denserks.a

#path to the L-BFGS-L solver
LBFGS= ../DENSERKS/examples/lbfgs.f

MODPATH= ../DENSERKS/lib/modules/

default: all

all: lib vanderpol

debug:
	$(MAKE) FFLAGS="-g -O0 -ffree-line-length-none"
	
lib:
	( cd ../DENSERKS/src; $(MAKE) )

vanderpol: vanderpol.f90 
	$(FC) $(FFLAGS) $< -o $@ $(DENSERKLIB) $(MODSWITCH) $(MODPATH)

clean:
	( cd ../DENSERKS/src; $(MAKE) clean )
	find . -maxdepth 1 -type f \
		! -name '*.f90' \
		! -name 'Makefile' \
		! -name 'makefile' \
		-exec rm -f {} +
	rm -f *.o *~
	
purge:
	( cd ../DENSERKS/src; $(MAKE) purge )
	
help:
	@echo "usage: make ?"

