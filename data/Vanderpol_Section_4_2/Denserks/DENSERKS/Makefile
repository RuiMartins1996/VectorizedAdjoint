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

#~~~> PGF90 (Portland Group Compiler)
#~~~> Optimize and generate code for the Pentium4 architecture
# export FC=pgf90
# export FFLAGS=-fast -tp=piv -Mcache_align -Mflushz -Minform=warn
# export MODSWITCH=-module

export FC=gfortran
#export FFLAGS= -Ofast -ffree-line-length-none $(LS_CONFIG)
export FFLAGS= -g -ffree-line-length-none $(LS_CONFIG)
export MODSWITCH=-I

#~~~> Lahey95 Compiler
#~~~> Optimize for SSE instruction set and Pentium4 architecture
#export FC=lf95
#export FFLAGS=--o2 --sse2 --tp4 --warn
#export MODSWITCH=-M

#~~~> Debugging flags for Lahey's lf95.
#~~~> Their use will result in a massive slowdown of the generated binaries (even 100-fold)
#     but they prove very useful for debugging purposes.
# export FFLAGS=-Cpp --chk a,e,s,u --pca --ap -O0 -g --trap --trace --chkglobal

#~~~> g95 Fortran Compiler
# export FC=g95
#~~~> Enabling optimizations with g95 will cause the Convection-Diffusion driver optimization
#     routine to fail (i.e., diverge); disable optimizations (remove -O switch) and recompile
#     to get the correct program behavior. 
# export FFLAGS=-O
# export MODSWITCH=-I

default: all

all:
	( cd ./src; $(MAKE) )
	( cd ./examples; $(MAKE) )

src:
	( cd ./src; $(MAKE) )

examples:
	( cd ./examples; $(MAKE) )

clean:
	( cd ./src; $(MAKE) clean )
	( cd ./examples; $(MAKE) clean )
	
purge:
	( cd ./src; $(MAKE) purge )
	( cd ./examples; $(MAKE) purge )
	
help:
	@$(ECHO) "usage: make ?"
