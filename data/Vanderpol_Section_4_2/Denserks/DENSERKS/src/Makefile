#Source makefile

#RANLIB=ranlib
#ARCH=ar cr

default: ../lib/denserks.a

#
FWD=erk.o
#
TLM=erk_tlm.o
#
ADJ=erk_adj.o erk_adjdr.o erk_adjdr_M.o
#
SOA=erk_soa.o erk_soadr.o erk_soadr_M.o
#
TAPES=erk_tapes.o erk_tlm_tapes.o
#
UTILS=erk_utils.o
#
MODULES=*.mod
#

erk_utils.o: erk_utils.f90
	$(FC) $(FFLAGS) -c $<

#-o $@

erk_tapes.o: erk_tapes.f90
	$(FC) $(FFLAGS) -c $<

#-o $@

erk_tlm_tapes.o: erk_tlm_tapes.f90
	$(FC) $(FFLAGS) -c $<
#-o $@

erk.o: erk.f90
	$(FC) $(FFLAGS) -c $<

#-o $@

erk_tlm.o: erk_tlm.f90
	$(FC) $(FFLAGS) -c $<

#-o $@

erk_adj.o: erk_adj.f90
	$(FC) $(FFLAGS) -c $<

# -o $@

erk_adjdr.o: erk_adjdr.f90
	$(FC) $(FFLAGS) -c $<

# -o $@

erk_adjdr_M.o: erk_adjdr_M.f90
	$(FC) $(FFLAGS) -c $<

# -o $@

erk_soa.o: erk_soa.f90
	$(FC) $(FFLAGS) -c $<

#-o $@

erk_soadr.o: erk_soadr.f90
	$(FC) $(FFLAGS) -c $<

#-o $@

erk_soadr_M.o: erk_soadr_M.f90
	$(FC) $(FFLAGS) -c $<

#-o $@

../lib/denserks.a: $(UTILS) $(TAPES) $(FWD) $(TLM) $(ADJ) $(SOA)
	$(ARCH) ../lib/denserks.a $^
	$(RANLIB) ../lib/denserks.a
	mv *.mod ../lib/modules/

purge: clean
	rm -f ../lib/denserks.a
	rm -f ../lib/modules/*.mod

clean:
	rm -f *.o *~
