#Makefile for the DENSERKS examples (driver programs)

default: dr_erk_soa_chk_optc \
	dr_erk_adj_chk_optc \
	dr_erk_adj_chk_arstf \
	dr_erk_adj_chk_arstf_M \
	dr_erk_adj_chk_vdp \
	dr_erk_adj_chk_cdiff \
	dr_erk_soa_chk_vdp \
	dr_erk_soa_chk_vdp_M \
	dr_erk_tlm_arstf \
	dr_erk_tlm_vdp \
	dr_erk_arstf \
	dr_erk_vdp

#path to the DENSERKS library
DENSERKLIB=../lib/denserks.a

#path to the L-BFGS-L solver
LBFGS=./lbfgs.f

MODPATH=../lib/modules/

dr_erk_soa_chk_optc: dr_erk_soa_chk_optc.f90
	$(FC) $(FFLAGS) $< -o $@ $(DENSERKLIB) $(MODSWITCH) $(MODPATH)
	mv dr_erk_soa_chk_optc ../bin

dr_erk_adj_chk_optc: dr_erk_adj_chk_optc.f90
	$(FC) $(FFLAGS) $< -o $@ $(DENSERKLIB) $(MODSWITCH) $(MODPATH)
	mv dr_erk_adj_chk_optc ../bin

dr_erk_adj_chk_arstf: dr_erk_adj_chk_arstf.f90 
	$(FC) $(FFLAGS) $< -o $@ $(DENSERKLIB) $(MODSWITCH) $(MODPATH)
	mv dr_erk_adj_chk_arstf ../bin

dr_erk_adj_chk_arstf_M: dr_erk_adj_chk_arstf_M.f90
	$(FC) $(FFLAGS) $< -o $@ $(DENSERKLIB) $(MODSWITCH) $(MODPATH)
	mv dr_erk_adj_chk_arstf_M ../bin

dr_erk_adj_chk_vdp: dr_erk_adj_chk_vdp.f90
	$(FC) $(FFLAGS) $< -o $@ $(DENSERKLIB) $(MODSWITCH) $(MODPATH)
	mv dr_erk_adj_chk_vdp ../bin

dr_erk_adj_chk_cdiff: dr_erk_adj_chk_cdiff.f90
	$(FC) $(FFLAGS) $< -o $@ $(DENSERKLIB) $(LBFGS) $(MODSWITCH) $(MODPATH)
	mv dr_erk_adj_chk_cdiff ../bin

dr_erk_soa_chk_vdp: dr_erk_soa_chk_vdp.f90
	$(FC) $(FFLAGS) $< -o $@ $(DENSERKLIB) $(MODSWITCH) $(MODPATH)
	mv dr_erk_soa_chk_vdp ../bin

dr_erk_soa_chk_vdp_M: dr_erk_soa_chk_vdp_M.f90
	$(FC) $(FFLAGS) $< -o $@ $(DENSERKLIB) $(MODSWITCH) $(MODPATH)
	mv dr_erk_soa_chk_vdp_M ../bin

dr_erk_tlm_arstf: dr_erk_tlm_arstf.f90
	$(FC) $(FFLAGS) $< -o $@ $(DENSERKLIB) $(MODSWITCH) $(MODPATH)
	mv dr_erk_tlm_arstf ../bin

dr_erk_tlm_vdp: dr_erk_tlm_vdp.f90
	$(FC) $(FFLAGS) $< -o $@ $(DENSERKLIB) $(MODSWITCH) $(MODPATH)
	mv dr_erk_tlm_vdp ../bin

dr_erk_vdp: dr_erk_vdp.f90
	$(FC) $(FFLAGS) $< -o $@ $(DENSERKLIB) $(MODSWITCH) $(MODPATH)
	mv dr_erk_vdp ../bin

dr_erk_arstf: dr_erk_arstf.f90
	$(FC) $(FFLAGS) $< -o $@ $(DENSERKLIB) $(MODSWITCH) $(MODPATH)
	mv dr_erk_arstf ../bin

purge: clean

clean:
	rm -f *.o *.mod *~
	rm -f ../bin/chk_T.dat \
	../bin/chk_Y.dat \
	../bin/chk_DY.dat \
	../bin/chk_H.dat \
	../bin/iterate.dat
	rm -f ../bin/dr_erk_soa_chk_optc \
	../bin/dr_erk_adj_chk_optc \
	../bin/dr_erk_adj_chk_arstf \
	../bin/dr_erk_adj_chk_arstf_M \
	../bin/dr_erk_adj_chk_vdp \
	../bin/dr_erk_adj_chk_cdiff \
	../bin/dr_erk_soa_chk_vdp \
	../bin/dr_erk_soa_chk_vdp_M \
	../bin/dr_erk_tlm_arstf \
	../bin/dr_erk_tlm_vdp \
	../bin/dr_erk_arstf \
	../bin/dr_erk_vdp
