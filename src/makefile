#                 ****************************************
#                               HydroBase 3 
#                 ****************************************
#                    Make file for HydroBase 3 Software
#
#                             Ruth Gorski Curry
#                      Woods Hole Oceanographic Institution
#                           Woods Hole, MA  02543
#
#
#-----------------------------------------------------------------------------
#  define these paths before making anything ...
#-----------------------------------------------------------------------------
INCLUDE= .
NETCDFLIB = /usr/local/netcdf_3/lib
NETCDFINC = /usr/local/netcdf_3/include
BINDIR = ../bin
#
#-----------------------------------------------------------------------------
DEBUG= -g   # set to -g to turn on debugging and set OPTIMIZE to a blank
OPTIMIZE =  #  optimization flags
LDEBUG=     #
#
# these define the compiler commands, flags and system libraries needed ...
#   appropriate options for some systems have been included ...
#   UNcomment the set of options appropriate to your machine
#
##
CC= gcc                       # invoke your C compiler
CFLAGS=-c -I$(INCLUDE)        # compile, include  
LINK= gcc                     # linker 
LDFLAGS =                     # any linker flags
LIBS=  -lm                 # math  libs
#LIBS= -lc  -lm                 # math, c libs
#
#-------------------------------------------------------------------------------
#       !! STOP EDITING HERE, THE REST IS FIXED !!
#-----------------------------------------------------------------------------
#
#
CDF_H = $(NETCDFINC)/netcdf.h
CDF_L = -L$(NETCDFLIB) -lnetcdf 
GRIDS_H = $(INCLUDE)/hb_grids.h
GAMMA_H = $(INCLUDE)/hb_gamma.h
HYDRO_H = $(INCLUDE)/hydrobase.h
HYDRO_CDF_H = $(INCLUDE)/hydro_cdf.h
MEMORY_H = $(INCLUDE)/hb_memory.h
NR_H =  $(INCLUDE)/hb_nr.h
OI3D_H = $(INCLUDE)/hb_oi3d.h
PATHS_H = $(INCLUDE)/hb_paths.h
SIGBINS_H = $(INCLUDE)/hb_sigbins.h
VGRAMS_H = $(INCLUDE)/hb_variograms.h
ZGRID_H = $(INCLUDE)/zgrid.h
#

PROG_SRC =  hb_bin3d.c hb_columns.c  hb_distdepth.c  hb_extract.c hb_random.c hb_ncfg3d.c  hb_foverh.c hb_getbottom.c hb_getpos.c   hb_grid2xyz.c hb_gridsection.c hb_gridsurf2d.c hb_layerav.c hb_msq10bounds.c  hb_mssort.c hb_mysort.c hb_nc2asc.c hb_ncinfo.c hb_ncmerge.c hb_ncoi3d.c hb_ncslice.c  hb_ncsmooth3d.c hb_ncsurf2d.c hb_propcalc.c hb_prseries.c hb_section.c hb_siftlevs.c hb_siftsta.c hb_smoothsurf2d.c  hb_stationsort.c  hb_surf2d.c hb_surfdiff2d.c hb_topo2xyz.c  hb_toposlice.c hb_updatehdr.c hb_xyprop.c hb_xyzprop.c woa_convert.c woce_btl2hb.c  woce_ctd2hb.c  wod_convert.c 
PROG_OBJ = ${PROG_SRC:.c=.o}
PROGRAMS = ${PROG_SRC:.c=}
UTILS_SRC = eos80.c hydro_cdf.c gamma_subs.c grid_subs.c hydro_utils.c julconv.c matrix_subs.c mem_subs.c oi3d_utils.c phyprops.c prop_subs.c topo_subs.c vgram_subs.c zgrid.c
UTILS_OBJ = ${UTILS_SRC:.c=.o}
OTHER_SRC = argo_convert.c  bats_ctdconvert.c  bats_btlconvert.c hb_add_ladcp.c hb_build_oi_parms.c hb_cond.c hb_coriolis.c hb_getdist.c hb_gettrack.c hb_ncfg3d.c hb_saltshift.c hb_theta.c  hb_tssig.c hb_vrotate.c  woce_btl_qualchk.c
OTHER_OBJ = ${OTHER_SRC:.c=.o}
OTHER = ${OTHER_SRC:.c=}
QC_SRC = hb_statchk_ts.c hb_statfit_ts.c hb_rangechk_ts.c hb_rangechk_to.c hb_saltshift.c  hb_statfit_to.c hb_statchk_to.c 
QC_OBJ = ${QC_SRC:.c=.o}
QC = ${QC_SRC:.c=}
#-----------------------------------------------------------------------------
#
hydrobase:	$(PROGRAMS)
#
qc: 		$(QC)
#
other:          $(OTHER)
#
install:	hydrobase 
	chmod 755 $(PROGRAMS)
	mv $(PROGRAMS) $(BINDIR)
#
install_qc:	qc 
	chmod 755 $(QC)
	mv $(QC) $(BINDIR)
#
install_other:	other 
	chmod 755 $(OTHER)
	mv $(OTHER) $(BINDIR)
#
clean:
	rm  -f $(PROG_OBJ) $(UTILS_OBJ) $(PROGRAMS) 
clean_all:
	rm  -f $(PROG_OBJ) $(UTILS_OBJ) $(PROGRAMS) $(BATH_OBJ) $(BATH) $(QC_OBJ) $(QC) $(OTHER_OBJ) $(OTHER)
#-----------------------------------------------------------------------------
#
#
#   PROGRAM DEPENDENCIES
#
hb_bin3d:	hb_bin3d.o hydro_cdf.o prop_subs.o hydro_utils.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_bin3d.o hydro_cdf.o prop_subs.o phyprops.o eos80.o hydro_utils.o   $(CDF_L) -o hb_bin3d $(LIBS) 
#
hb_columns: 	hb_columns.o phyprops.o prop_subs.o hydro_utils.o gamma_subs.o eos80.o
	$(LINK) $(LDFLAGS)  hb_columns.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o $(CDF_L) -o hb_columns $(LIBS)
#
hb_coriolis: 	hb_coriolis.o phyprops.o  eos80.o
	$(LINK) $(LDFLAGS)  hb_coriolis.o phyprops.o eos80.o  -o hb_coriolis $(LIBS)
#
hb_distdepth:  	hb_distdepth.o topo_subs.o grid_subs.o mem_subs.o
	$(LINK) $(LDFLAGS) hb_distdepth.o topo_subs.o grid_subs.o mem_subs.o $(CDF_L) -o hb_distdepth $(LIBS)
#
hb_extract:	hb_extract.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_extract.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_extract $(LIBS)
#
hb_random:	hb_random.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_random.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_random $(LIBS)
#
hb_findblanks:		hb_findblanks.o 
	$(LINK) $(LDFLAGS) hb_findblanks.o $(LIBS) -o hb_findblanks 
#
hb_foverh:	hb_foverh.o mem_subs.o phyprops.o eos80.o 
	$(LINK) $(LDFLAGS) hb_foverh.o  phyprops.o eos80.o  mem_subs.o $(LIBS) -o hb_foverh 
#
hb_getbottom: 	hb_getbottom.o phyprops.o prop_subs.o hydro_utils.o gamma_subs.o  eos80.o
	$(LINK) $(LDFLAGS)  hb_getbottom.o phyprops.o prop_subs.o hydro_utils.o eos80.o  gamma_subs.o  $(CDF_L)  -o hb_getbottom $(LIBS)
#
hb_getdist:	hb_getdist.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_getdist.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_getdist $(LIBS)
#
hb_getpos:	hb_getpos.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_getpos.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_getpos $(LIBS)
#
hb_gettrack:  	hb_gettrack.o 
	$(LINK) $(LDFLAGS) hb_gettrack.o  -o hb_gettrack $(LIBS)
#
hb_grid2xyz:	hb_grid2xyz.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_grid2xyz.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_grid2xyz $(LIBS)
#
hb_gridsection:	hb_gridsection.o grid_subs.o mem_subs.o
	$(LINK) $(LDFLAGS) hb_gridsection.o grid_subs.o  mem_subs.o $(LIBS) -o hb_gridsection 
#
hb_gridsurf2d: 	hb_gridsurf2d.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o
	$(LINK) $(LDFLAGS)  hb_gridsurf2d.o phyprops.o prop_subs.o  hydro_utils.o gamma_subs.o  eos80.o  $(CDF_L) -o hb_gridsurf2d $(LIBS)
#
hb_layerav: 	hb_layerav.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o
	$(LINK) $(LDFLAGS)  hb_layerav.o phyprops.o prop_subs.o  hydro_utils.o eos80.o gamma_subs.o $(CDF_L) $(LIBS)  -o hb_layerav
#
hb_layerinteg: 	hb_layerinteg.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o
	$(LINK) $(LDFLAGS)  hb_layerinteg.o phyprops.o prop_subs.o  hydro_utils.o eos80.o gamma_subs.o $(CDF_L) $(LIBS)  -o hb_layerinteg
#
hb_min_max:		hb_min_max.o 
	$(LINK) $(LDFLAGS) hb_min_max.o $(LIBS) -o hb_min_max
hb_msq10bounds:	hb_msq10bounds.o
	$(LINK) $(LDFLAGS) hb_msq10bounds.o -o hb_msq10bounds $(LIBS)
#
hb_mssort:	hb_mssort.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_mssort.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_mssort $(LIBS)
#
hb_mysort:	hb_mysort.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_mysort.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_mysort $(LIBS)
#
hb_nc2asc:	hb_nc2asc.o hydro_utils.o hydro_cdf.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_nc2asc.o hydro_utils.o hydro_cdf.o prop_subs.o phyprops.o  eos80.o $(CDF_L) -o hb_nc2asc $(LIBS)
#
hb_ncinfo:	hb_ncinfo.o hydro_cdf.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_ncinfo.o hydro_cdf.o prop_subs.o phyprops.o  eos80.o $(CDF_L) -o hb_ncinfo $(LIBS)
#
hb_ncmerge:	hb_ncmerge.o hydro_utils.o hydro_cdf.o prop_subs.o phyprops.o eos80.o  mem_subs.o grid_subs.o topo_subs.o
	$(LINK) $(LDFLAGS) hb_ncmerge.o hydro_utils.o hydro_cdf.o prop_subs.o phyprops.o  eos80.o  mem_subs.o $(CDF_L) grid_subs.o topo_subs.o -o hb_ncmerge $(LIBS)
#
hb_ncoi3d:	hb_ncoi3d.o hydro_utils.o hydro_cdf.o prop_subs.o phyprops.o eos80.o matrix_subs.o mem_subs.o oi3d_utils.o grid_subs.o topo_subs.o vgram_subs.o
	$(LINK) $(LDFLAGS)  hb_ncoi3d.o mem_subs.o phyprops.o prop_subs.o hydro_cdf.o hydro_utils.o eos80.o  grid_subs.o topo_subs.o vgram_subs.o oi3d_utils.o matrix_subs.o $(CDF_L)  -o hb_ncoi3d $(LIBS)
#
hb_ncslice: 	hb_ncslice.o mem_subs.o phyprops.o prop_subs.o hydro_cdf.o hydro_utils.o eos80.o gamma_subs.o grid_subs.o topo_subs.o
	$(LINK) $(LDFLAGS)  hb_ncslice.o mem_subs.o phyprops.o prop_subs.o hydro_cdf.o hydro_utils.o eos80.o  gamma_subs.o grid_subs.o topo_subs.o $(CDF_L)  -o hb_ncslice $(LIBS)
#
hb_ncsmooth3d:	hb_ncsmooth3d.o grid_subs.o hydro_utils.o hydro_cdf.o prop_subs.o phyprops.o eos80.o mem_subs.o
	$(LINK) $(LDFLAGS) hb_ncsmooth3d.o grid_subs.o  hydro_utils.o hydro_cdf.o  prop_subs.o phyprops.o eos80.o mem_subs.o $(CDF_L)  -o hb_ncsmooth3d $(LIBS)
#
hb_ncsurf2d: 	hb_ncsurf2d.o phyprops.o prop_subs.o hydro_cdf.o hydro_utils.o eos80.o gamma_subs.o
	$(LINK) $(LDFLAGS)  hb_ncsurf2d.o phyprops.o prop_subs.o hydro_cdf.o hydro_utils.o gamma_subs.o  eos80.o   $(CDF_L)  -o hb_ncsurf2d $(LIBS)
#
hb_propcalc: 	hb_propcalc.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o 
	$(LINK) $(LDFLAGS)  hb_propcalc.o phyprops.o prop_subs.o hydro_utils.o eos80.o  gamma_subs.o $(CDF_L) -o hb_propcalc $(LIBS)
#
hb_prseries: 	hb_prseries.o phyprops.o prop_subs.o hydro_utils.o eos80.o
	$(LINK) $(LDFLAGS)  hb_prseries.o phyprops.o prop_subs.o hydro_utils.o eos80.o   $(CDF_L)  -o hb_prseries $(LIBS)
#
hb_saltshift: 	hb_saltshift.o phyprops.o prop_subs.o hydro_utils.o eos80.o
	$(LINK) $(LDFLAGS)  hb_saltshift.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o hb_saltshift $(LIBS)
#
hb_section: 	hb_section.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o
	$(LINK) $(LDFLAGS)  hb_section.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o  $(CDF_L)  -o hb_section $(LIBS)
#
hb_siftlevs:	hb_siftlevs.o hydro_utils.o prop_subs.o phyprops.o eos80.o gamma_subs.o 
	$(LINK) $(LDFLAGS) hb_siftlevs.o hydro_utils.o prop_subs.o phyprops.o eos80.o  gamma_subs.o  $(CDF_L) -o hb_siftlevs $(LIBS)
#
hb_siftsta:	hb_siftsta.o hydro_utils.o prop_subs.o phyprops.o eos80.o gamma_subs.o
	$(LINK) $(LDFLAGS) hb_siftsta.o hydro_utils.o prop_subs.o phyprops.o eos80.o gamma_subs.o  $(CDF_L) -o hb_siftsta $(LIBS)
#
hb_sigma:	hb_sigma.o  phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_sigma.o phyprops.o eos80.o -o hb_sigma $(LIBS)
#
hb_sigma_eval:	hb_sigma_eval.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_sigma_eval.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_sigma_eval $(LIBS)
#
hb_sigma024_eval:	hb_sigma024_eval.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_sigma024_eval.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_sigma024_eval $(LIBS)
#
hb_smoothsurf2d:	hb_smoothsurf2d.o  grid_subs.o mem_subs.o
	$(LINK) $(LDFLAGS) hb_smoothsurf2d.o grid_subs.o mem_subs.o -o hb_smoothsurf2d $(LIBS)
#
hb_stationsort:	hb_stationsort.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_stationsort.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_stationsort $(LIBS)
#
hb_stats1d:	hb_stats1d.o 
	$(LINK) $(LDFLAGS)  hb_stats1d.o  -o hb_stats1d $(LIBS)
#
hb_surf2d: 	hb_surf2d.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o
	$(LINK) $(LDFLAGS)  hb_surf2d.o phyprops.o prop_subs.o  hydro_utils.o eos80.o gamma_subs.o $(CDF_L) $(LIBS)  -o hb_surf2d
#
hb_surfdiff2d:	hb_surfdiff2d.o 
	$(LINK) $(LDFLAGS) hb_surfdiff2d.o $(LIBS) -o hb_surfdiff2d 
#
hb_theta:	hb_theta.o  phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_theta.o phyprops.o eos80.o -o hb_theta $(LIBS)
#
hb_topo2xyz:	hb_topo2xyz.o topo_subs.o grid_subs.o mem_subs.o
	$(LINK) $(LDFLAGS) hb_topo2xyz.o topo_subs.o grid_subs.o mem_subs.o $(CDF_L) $(LIBS) -o hb_topo2xyz 
#
hb_toposlice:	hb_toposlice.o mem_subs.o  phyprops.o prop_subs.o hydro_utils.o eos80.o topo_subs.o grid_subs.o
	$(LINK) $(LDFLAGS) hb_toposlice.o  mem_subs.o  phyprops.o prop_subs.o hydro_utils.o eos80.o topo_subs.o grid_subs.o $(CDF_L) $(LIBS) -o hb_toposlice 
#
hb_tssig: 	hb_tssig.o phyprops.o prop_subs.o hydro_utils.o eos80.o
	$(LINK) $(LDFLAGS)  hb_tssig.o phyprops.o prop_subs.o  hydro_utils.o eos80.o  $(CDF_L) $(LIBS)  -o hb_tssig
#
hb_updatehdr:	hb_updatehdr.o phyprops.o prop_subs.o hydro_utils.o eos80.o topo_subs.o grid_subs.o mem_subs.o
	$(LINK) $(LDFLAGS)  hb_updatehdr.o phyprops.o prop_subs.o  hydro_utils.o eos80.o topo_subs.o grid_subs.o mem_subs.o $(CDF_L)   -o hb_updatehdr $(LIBS)
#
hb_volumetric: 	hb_volumetric.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o
	$(LINK) $(LDFLAGS)  hb_volumetric.o phyprops.o prop_subs.o  hydro_utils.o eos80.o gamma_subs.o $(CDF_L) $(LIBS)  -o hb_volumetric
#
hb_volumetric_anom: 	hb_volumetric_anom.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o
	$(LINK) $(LDFLAGS)  hb_volumetric_anom.o phyprops.o prop_subs.o  hydro_utils.o eos80.o gamma_subs.o $(CDF_L) $(LIBS)  -o hb_volumetric_anom
#
hb_xyprop:	hb_xyprop.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o 
	$(LINK) $(LDFLAGS)  hb_xyprop.o phyprops.o prop_subs.o  hydro_utils.o eos80.o gamma_subs.o $(CDF_L)  -o hb_xyprop $(LIBS)
#
hb_xyzprop:	hb_xyzprop.o phyprops.o prop_subs.o hydro_utils.o eos80.o gamma_subs.o 
	$(LINK) $(LDFLAGS)  hb_xyzprop.o phyprops.o prop_subs.o  hydro_utils.o eos80.o gamma_subs.o $(CDF_L)  -o hb_xyzprop $(LIBS)
#
ices_convert:	ices_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o filters.o mem_subs.o
	$(LINK) $(LDFLAGS)  ices_convert.o  phyprops.o prop_subs.o hydro_utils.o eos80.o filters.o mem_subs.o $(CDF_L) -o ices_convert $(LIBS)
#
sac_convert:	sac_convert.o hydro_cdf.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  sac_convert.o  hydro_cdf.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o sac_convert $(LIBS)
#
seabird_btl_convert:	seabird_btl_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  seabird_btl_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o seabird_btl_convert $(LIBS)
#
seabird_dcc_convert:	seabird_dcc_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  seabird_dcc_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o seabird_dcc_convert $(LIBS)
#
woa_convert:	woa_convert.o hydro_cdf.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) woa_convert.o hydro_cdf.o hydro_utils.o prop_subs.o phyprops.o eos80.o $(CDF_L) -o woa_convert $(LIBS)
#
woce_btl2hb:	woce_btl2hb.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  woce_btl2hb.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o woce_btl2hb $(LIBS)
#
woce_ctd2hb:	woce_ctd2hb.o hydro_utils.o prop_subs.o phyprops.o eos80.o mem_subs.o
	$(LINK) $(LDFLAGS)  woce_ctd2hb.o phyprops.o prop_subs.o hydro_utils.o eos80.o  mem_subs.o  $(CDF_L) -o woce_ctd2hb $(LIBS)
#
wod_convert:	wod_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  wod_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o  -o wod_convert $(LIBS)
#
#
#
#
#   OBJECT DEPENDENCIES
#
eos80.o:	eos80.c $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) eos80.c
#
grid_subs.o:	grid_subs.c $(GRIDS_H) $(MEMORY_H)
	$(CC)  $(CFLAGS) $(DEBUG) grid_subs.c
#
gamma_subs.o:	gamma_subs.c $(GAMMA_H) $(HYDRO_H) $(CDF_H)
	$(CC) $(CFLAGS)  -I$(NETCDFINC) $(DEBUG) gamma_subs.c
#
hydro_cdf.o:	hydro_cdf.c  $(HYDRO_H) $(HYDRO_CDF_H) $(CDF_H)
	$(CC) $(CFLAGS)  -I$(NETCDFINC) $(DEBUG) hydro_cdf.c
#
hydro_utils.o:	hydro_utils.c $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) hydro_utils.c
#
matrix_subs.o:	matrix_subs.c $(NR_H)
	$(CC) $(CFLAGS) $(DEBUG) matrix_subs.c
#
mem_subs.o:	mem_subs.c $(MEMORY_H)
	$(CC) $(CFLAGS) $(DEBUG) mem_subs.c
#
oi3d_utils.o:	oi3d_utils.c $(OI3D_H) $(CDF_H)
	$(CC)  $(CFLAGS) -I$(NETCDFINC)  $(DEBUG) oi3d_utils.c
#
phyprops.o:	phyprops.c $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) phyprops.c
#
prop_subs.o:	prop_subs.c $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) prop_subs.c
#
topo_subs.o:	topo_subs.c $(GRIDS_H)  $(CDF_H)
	$(CC)  $(CFLAGS) -I$(NETCDFINC) $(DEBUG) topo_subs.c
#
vgram_subs.o:   vgram_subs.c $(VGRAMS_H) $(GRIDS.H) $(MEMORY_H) 
	$(CC) $(CFLAGS) $(DEBUG) vgram_subs.c
# 
zgrid.o:	zgrid.c $(ZGRID_H) $(MEMORY_H)
	$(CC)  $(CFLAGS) $(DEBUG) zgrid.c
#
#################
hb_bin3d.o:	hb_bin3d.c  $(HYDRO_H) $(HYDRO_CDF_H) $(PATHS_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_bin3d.c
#
hb_columns.o:	hb_columns.c $(HYDRO_H) $(PATHS_H) $(GAMMA_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_columns.c
#
hb_coriolis.o:       hb_coriolis.c $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_coriolis.c
#
hb_distdepth.o:	hb_distdepth.c $(MEMORY_H) $(PATHS_H) $(GRIDS_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_distdepth.c
#
hb_extract.o:	hb_extract.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_extract.c
#
hb_random.o:	hb_random.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_random.c
#
hb_findblanks.o:	hb_findblanks.c 
	$(CC)  $(CFLAGS) $(DEBUG) hb_findblanks.c
#
hb_foverh.o:	hb_foverh.c  $(MEMORY_H) $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_foverh.c
#
hb_getbottom.o:	hb_getbottom.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_getbottom.c
#
hb_getdist.o:	hb_getdist.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_getdist.c
#
hb_getpos.o:	hb_getpos.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_getpos.c
#
hb_gettrack.o:	hb_gettrack.c $(GRIDS_H) 
	$(CC)  $(CFLAGS) $(DEBUG) hb_gettrack.c
#
hb_grid2xyz.o:	hb_grid2xyz.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_grid2xyz.c
#
hb_gridsection.o:	hb_gridsection.c  $(HYDRO_H) $(GRIDS_H) $(MEMORY_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_gridsection.c
#
hb_gridsurf2d.o:	hb_gridsurf2d.c   $(HYDRO_H) $(PATHS_H) $(GAMMA_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_gridsurf2d.c
#
hb_interpolate.o:	hb_interpolate.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_interpolate.c
#
hb_layerav.o:	hb_layerav.c   $(HYDRO_H) $(PATHS_H) $(GAMMA_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_layerav.c
#
hb_layerinteg.o:	hb_layerinteg.c   $(HYDRO_H) $(PATHS_H) $(GAMMA_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_layerinteg.c
#
hb_min_max.o:	hb_min_max.c 
	$(CC)  $(CFLAGS) $(DEBUG) hb_min_max.c
#
hb_msq10bounds.o:	hb_msq10bounds.c
	$(CC)  $(CFLAGS) $(DEBUG) hb_msq10bounds.c
#
hb_mssort.o:	hb_mssort.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_mssort.c
#
hb_mysort.o:	hb_mysort.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_mysort.c
#
hb_nc2asc.o:       hb_nc2asc.c $(HYDRO_CDF_H) $(HYDRO_H) $(CDF_H) 
	$(CC) $(CFLAGS) $(DEBUG)  -I$(NETCDFINC)  hb_nc2asc.c
#
hb_ncinfo.o:       hb_ncinfo.c $(HYDRO_CDF_H) $(CDF_H) 
	$(CC) $(CFLAGS) $(DEBUG)  -I$(NETCDFINC)  hb_ncinfo.c
#
hb_ncmerge.o:       hb_ncmerge.c $(HYDRO_CDF_H) $(HYDRO_H) $(MEMORY_H) $(PATHS_H)  $(GRIDS_H) $(CDF_H)
	$(CC) $(CFLAGS) $(DEBUG)  -I$(NETCDFINC)  hb_ncmerge.c
#
hb_ncoi3d.o:	hb_ncoi3d.c $(CDF_H)  $(HYDRO_H) $(HYDRO_CDF_H) $(MEMORY_H) $(NC_H) $(PATHS_H) $(GRIDS_H) $(VGRAMS_H) $(OI3D_H)
	$(CC) $(CFLAGS) $(DEBUG)  -I$(NETCDFINC)  hb_ncoi3d.c
#
hb_ncslice.o:	hb_ncslice.c  $(CDF_H)  $(HYDRO_H) $(HYDRO_CDF_H) $(MEMORY_H) $(PATHS_H) $(GAMMA_H)   $(GRIDS_H)
	$(CC) $(CFLAGS) $(DEBUG)  -I$(NETCDFINC)  hb_ncslice.c
#
hb_ncsmooth3d.o:	hb_ncsmooth3d.c $(HYDRO_H) $(HYDRO_CDF_H) $(GRIDS_H) $(MEMORY_H) $(CDF_H)
	$(CC)  $(CFLAGS) $(DEBUG)   -I$(NETCDFINC) hb_ncsmooth3d.c
#
hb_ncsurf2d.o:	hb_ncsurf2d.c   $(HYDRO_H) $(HYDRO_CDF_H) $(PATHS_H) $(GAMMA_H)  $(CDF_H)
	$(CC) $(CFLAGS) $(DEBUG)  -I$(NETCDFINC)  hb_ncsurf2d.c
#
hb_prseries.o:	hb_prseries.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_prseries.c
#
hb_propcalc.o:	hb_propcalc.c $(HYDRO_H) $(PATHS_H) $(GAMMA_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_propcalc.c
#
hb_saltshift.o:	hb_saltshift.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_saltshift.c
#
hb_section.o:	hb_section.c $(HYDRO_H) $(PATHS_H) $(GAMMA_H)  
	$(CC)  $(CFLAGS) $(DEBUG) hb_section.c
#
hb_siftlevs.o:	hb_siftlevs.c  $(HYDRO_H) $(PATHS_H) $(GAMMA_H) 
	$(CC) $(CFLAGS) $(DEBUG) hb_siftlevs.c
#
hb_siftsta.o:	hb_siftsta.c  $(HYDRO_H) $(PATHS_H) $(GAMMA_H) 
	$(CC) $(CFLAGS) $(DEBUG) hb_siftsta.c
#
hb_sigma.o:	hb_sigma.c  $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_sigma.c
#
hb_sigma_eval.o:	hb_sigma_eval.c  $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_sigma_eval.c
#
hb_sigma024_eval.o:	hb_sigma024_eval.c  $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_sigma024_eval.c
#
hb_smoothsurf2d.o:	hb_smoothsurf2d.c $(MEMORY_H) $(GRIDS_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_smoothsurf2d.c
#
hb_stationsort.o:	hb_stationsort.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_stationsort.c
#
hb_stats1d.o:	hb_stats1d.c   
	$(CC) $(CFLAGS) $(DEBUG) hb_stats1d.c
#
hb_surf2d.o:	hb_surf2d.c   $(HYDRO_H) $(PATHS_H) $(GAMMA_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_surf2d.c
#
hb_surfdiff2d.o:	hb_surfdiff2d.c 
	$(CC)  $(CFLAGS) $(DEBUG) hb_surfdiff2d.c
#
hb_theta.o:	hb_theta.c  $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_theta.c
#
hb_topo2xyz.o:	hb_topo2xyz.c  $(MEMORY_H) $(CDF_H) $(GRIDS_H) $(PATHS_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_topo2xyz.c
#
hb_toposlice.o:	hb_toposlice.c  $(MEMORY_H)  $(HYDRO_H) $(PATHS_H) $(GRIDS_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_toposlice.c
#
hb_tssig.o:	hb_tssig.c   
	$(CC) $(CFLAGS) $(DEBUG) hb_tssig.c
#
hb_updatehdr.o:	hb_updatehdr.c $(HYDRO_H) $(CDF_H) $(GRIDS_H) $(PATHS_H)
	$(CC) $(CFLAGS) -I$(NETCDFINC)  $(DEBUG) hb_updatehdr.c
#
hb_volumetric.o:	hb_volumetric.c   $(HYDRO_H) $(PATHS_H) $(GAMMA_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_volumetric.c
#
hb_volumetric_anom.o:	hb_volumetric_anom.c   $(HYDRO_H) $(PATHS_H) $(GAMMA_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_volumetric_anom.c
#
hb_xyprop.o:	hb_xyprop.c   $(HYDRO_H) $(PATHS_H) $(GAMMA_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_xyprop.c
#
hb_xyzprop.o:	hb_xyzprop.c   $(HYDRO_H) $(PATHS_H) $(GAMMA_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_xyzprop.c
#
woce_btl2hb.o:	woce_btl2hb.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) woce_btl2hb.c
#
woce_ctd2hb.o:	woce_ctd2hb.c $(HYDRO_H) 
	$(CC)  $(CFLAGS) $(DEBUG) woce_ctd2hb.c
#
woa_convert.o:	woa_convert.c $(HYDRO_H) $(HYDRO_CDF_H)
	$(CC)  $(CFLAGS) $(DEBUG) woa_convert.c
#
wod_convert.o:	wod_convert.c $(HYDRO_H) $(PATHS_H)
	$(CC)  $(CFLAGS) $(DEBUG) wod_convert.c
#
#
#
#
#
#
#    OTHER PROGRAMS
#
argo_convert:	argo_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o  julconv.o
	$(LINK) $(LDFLAGS)  argo_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o  julconv.o  -o argo_convert $(CDF_L) $(LIBS)
#
awi_convert:	awi_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o  mem_subs.o
	$(LINK) $(LDFLAGS)  awi_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o mem_subs.o  -o awi_convert $(LIBS)
#
bbsr_ctdconvert:     bbsr_ctdconvert.o hydro_utils.o phyprops.o eos80.o prop_subs.o 
	$(LINK) $(LDFLAGS) bbsr_ctdconvert.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o bbsr_ctdconvert  $(LIBS)
#
bats_ctdconvert:     bats_ctdconvert.o hydro_utils.o phyprops.o eos80.o prop_subs.o 
	$(LINK) $(LDFLAGS) bats_ctdconvert.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o bats_ctdconvert  $(LIBS)
#
bats_btlconvert:     bats_btlconvert.o hydro_utils.o phyprops.o eos80.o prop_subs.o 
	$(LINK) $(LDFLAGS) bats_btlconvert.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o bats_btlconvert  $(LIBS)
#
odvconvert:     odvconvert.o hydro_utils.o phyprops.o eos80.o prop_subs.o 
	$(LINK) $(LDFLAGS) odvconvert.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o odvconvert  $(LIBS)
#
med_convert:     med_convert.o hydro_utils.o phyprops.o eos80.o prop_subs.o 
	$(LINK) $(LDFLAGS) med_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o med_convert  $(LIBS)
#
coastal_sta:	coastal_sta.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) coastal_sta.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o coastal_sta $(LIBS)
#
csv_btlconvert:     csv_btlconvert.o hydro_utils.o phyprops.o eos80.o prop_subs.o mem_subs.o
	$(LINK) $(LDFLAGS) csv_btlconvert.o hydro_utils.o prop_subs.o phyprops.o eos80.o mem_subs.o -o csv_btlconvert  $(LIBS)
#
ctda_convert:	ctda_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o filters.o mem_subs.o
	$(LINK) $(LDFLAGS)  ctda_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o filters.o mem_subs.o  -o ctda_convert $(LIBS)
#
hb_build_oi_parms:	hb_build_oi_parms.o  hydro_utils.o prop_subs.o phyprops.o eos80.o oi3d_utils.o topo_subs.o grid_subs.o vgram_subs.o hydro_cdf.o mem_subs.o
	$(LINK) $(LDFLAGS) hb_build_oi_parms.o  hydro_utils.o prop_subs.o phyprops.o eos80.o oi3d_utils.o topo_subs.o grid_subs.o vgram_subs.o hydro_cdf.o  mem_subs.o -o hb_build_oi_parms $(CDF_L) $(LIBS)
#
hb_add_ladcp:	hb_add_ladcp.o  hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_add_ladcp.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_add_ladcp  $(LIBS)
#
hb_cond:	hb_cond.o  phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_cond.o phyprops.o eos80.o -o hb_cond $(LIBS)
#
hb_ncfg3d:	hb_ncfg3d.o zgrid.o hydro_cdf.o prop_subs.o hydro_utils.o phyprops.o eos80.o mem_subs.o topo_subs.o grid_subs.o
	$(LINK) $(LDFLAGS) hb_ncfg3d.o hydro_cdf.o prop_subs.o phyprops.o eos80.o hydro_utils.o grid_subs.o  zgrid.o  mem_subs.o topo_subs.o $(CDF_L) $(LIBS) -o hb_ncfg3d 
#
hb_get_month_of_smin:	hb_get_month_of_smin.o 
	$(LINK) $(LDFLAGS) hb_get_month_of_smin.o $(LIBS) -o hb_get_month_of_smin
#
hb_linew_avg_section:	hb_linew_avg_section.o hydro_utils.o prop_subs.o phyprops.o eos80.o  gamma_subs.o 
	$(LINK) $(LDFLAGS) hb_linew_avg_section.o  hydro_utils.o prop_subs.o  phyprops.o eos80.o gamma_subs.o $(CDF_L) -o hb_linew_avg_section $(LIBS)
#
hb_surfmerge:	hb_surfmerge.o 
	$(LINK) $(LDFLAGS) hb_surfmerge.o $(LIBS) -o hb_surfmerge 
#
hb_test_oi_parms:	hb_test_oi_parms.o  hydro_utils.o prop_subs.o phyprops.o eos80.o topo_subs.o grid_subs.o vgram_subs.o hydro_cdf.o mem_subs.o
	$(LINK) $(LDFLAGS) hb_test_oi_parms.o  hydro_utils.o prop_subs.o phyprops.o eos80.o topo_subs.o grid_subs.o vgram_subs.o hydro_cdf.o  mem_subs.o -o hb_test_oi_parms $(CDF_L) $(LIBS)
#
hb_vrotate:	hb_vrotate.o  hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_vrotate.o  hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_vrotate  $(LIBS)
#
ifremer_convert:     ifremer_convert.o hydro_utils.o phyprops.o eos80.o prop_subs.o julconv.o
	$(LINK) $(LDFLAGS) ifremer_convert.o hydro_utils.o prop_subs.o phyprops.o eos80.o  julconv.o $(CDF_L) -o ifremer_convert  $(LIBS)
#
igor_ctdconvert:     igor_ctdconvert.o hydro_utils.o phyprops.o eos80.o prop_subs.o mem_subs.o
	$(LINK) $(LDFLAGS) igor_ctdconvert.o hydro_utils.o prop_subs.o phyprops.o eos80.o mem_subs.o -o igor_ctdconvert  $(LIBS)
#
gamma_ex:	gamma_ex.o gamma_subs.o phyprops.o eos80.o 
	$(LINK) $(LDFLAGS) gamma_ex.o gamma_subs.o phyprops.o eos80.o  $(CDF_L) -o gamma_ex $(LIBS)
#
kiel_ctdconvert:	kiel_ctdconvert.o  hydro_utils.o prop_subs.o phyprops.o eos80.o 
	$(LINK) $(LDFLAGS)  kiel_ctdconvert.o phyprops.o prop_subs.o hydro_utils.o eos80.o -o kiel_ctdconvert $(LIBS)
#
sbe_ctdconvert:	sbe_ctdconvert.o  hydro_utils.o prop_subs.o phyprops.o eos80.o 
	$(LINK) $(LDFLAGS)  sbe_ctdconvert.o phyprops.o prop_subs.o hydro_utils.o eos80.o -o sbe_ctdconvert $(LIBS)
#
mmp_convert:	mmp_convert.o  hydro_utils.o prop_subs.o phyprops.o eos80.o mem_subs.o
	$(LINK) $(LDFLAGS)  mmp_convert.o phyprops.o prop_subs.o hydro_utils.o eos80.o mem_subs.o  -o mmp_convert $(LIBS)
#
odf_ctdconvert:	odf_ctdconvert.o hydro_utils.o prop_subs.o phyprops.o eos80.o filters.o mem_subs.o
	$(LINK) $(LDFLAGS)  odf_ctdconvert.o phyprops.o prop_subs.o hydro_utils.o eos80.o filters.o mem_subs.o  -o odf_ctdconvert $(LIBS)
#
pe_transport:	pe_transport.o  phyprops.o eos80.o  hydro_utils.o prop_subs.o 
	$(LINK) $(LDFLAGS) pe_transport.o hydro_utils.o prop_subs.o phyprops.o eos80.o   $(LIBS)  -o pe_transport 
#
soc_ctdconvert:	soc_ctdconvert.o hydro_utils.o prop_subs.o phyprops.o eos80.o  mem_subs.o
	$(LINK) $(LDFLAGS)  soc_ctdconvert.o phyprops.o prop_subs.o hydro_utils.o eos80.o  mem_subs.o  -o soc_ctdconvert $(LIBS)
yearsmooth:	yearsmooth.o
	$(CC) $(DEBUG) yearsmooth.o -o yearsmooth -lm
nc_add_att:     nc_add_att.o 
	$(LINK) $(LDFLAGS) nc_add_att.o $(CDF_L) -o nc_add_att $(LIBS)
#
#    OTHER OBJECT DEPENDENCIES
#
argo_convert.o:	argo_convert.c $(HYDRO_H) 
	$(CC) $(CFLAGS)  -I$(NETCDFINC)  $(DEBUG) argo_convert.c
#
awi_convert.o:	awi_convert.c $(HYDRO_H) $(MEMORY_H)
	$(CC) $(CFLAGS) $(DEBUG) awi_convert.c
#
bats_ctdconvert.o:	bats_ctdconvert.c $(HYDRO_H) 
	$(CC) $(CFLAGS) $(DEBUG) bats_ctdconvert.c
#
bats_btlconvert.o:	bats_btlconvert.c $(HYDRO_H) 
	$(CC) $(CFLAGS) $(DEBUG) bats_btlconvert.c
#
odvconvert.o:	odvconvert.c $(HYDRO_H) 
	$(CC) $(CFLAGS) $(DEBUG) odvconvert.c
#
med_convert.o:	med_convert.c $(HYDRO_H) 
	$(CC) $(CFLAGS) $(DEBUG) med_convert.c
#
nc_add_att.o:	nc_add_att.c  
	$(CC) $(CFLAGS) $(DEBUG) -I$(NETCDFINC) nc_add_att.c
#
hb_build_oi_parms.o:	hb_build_oi_parms.c $(HYDRO_H) $(CDF_H) $(HYDRO_CDF_H) $(MEMORY_H)  $(PATHS_H) $(GRIDS_H) $(OI3D_H) $(VGRAMS_H) 
	$(CC) $(CFLAGS)  -I$(NETCDFINC)  $(DEBUG) hb_build_oi_parms.c
#
hb_cond.o:	hb_cond.c  $(HYDRO_H)
	$(CC) $(CFLAGS) $(DEBUG) hb_cond.c
#
hb_ncfg3d.o:	hb_ncfg3d.c $(ZGRID_H) $(HYDRO_H) $(HYDRO_CDF_H) $(MEMORY_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_ncfg3d.c
#
hb_add_ladcp.o:	 hb_add_ladcp.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_add_ladcp.c
#
hb_get_month_of_smin.o:	hb_get_month_of_smin.c 
	$(CC)  $(CFLAGS) $(DEBUG) hb_get_month_of_smin.c
#
hb_linew_avg_section.o:	hb_linew_avg_section.c $(HYDRO_H)  $(GAMMA_H) $(PATHS_H) 
	$(CC)  $(CFLAGS) $(DEBUG) hb_linew_avg_section.c
#
hb_surfmerge.o:	hb_surfmerge.c 
	$(CC)  $(CFLAGS) $(DEBUG) hb_surfmerge.c
#
hb_test_oi_parms.o:	hb_test_oi_parms.c $(HYDRO_H) $(CDF_H) $(HYDRO_CDF_H) $(MEMORY_H)  $(PATHS_H) $(GRIDS_H) $(VGRAMS_H) 
	$(CC) $(CFLAGS)  -I$(NETCDFINC)  $(DEBUG) hb_test_oi_parms.c
#
hb_vrotate.o:	 hb_vrotate.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_vrotate.c
#
ifremer_convert.o:	ifremer_convert.c $(HYDRO_H) $(CDF_H) 
	$(CC) $(CFLAGS) $(DEBUG) -I$(NETCDFINC) ifremer_convert.c
#
igor_ctdconvert.o:	igor_ctdconvert.c $(HYDRO_H) $(FILTERS_H)
	$(CC) $(CFLAGS) $(DEBUG) igor_ctdconvert.c
#
gamma_ex.o:	gamma_ex.c $(GAMMA_H) $(PATHS_H)
	$(CC)  $(CFLAGS) $(DEBUG) gamma_ex.c
#
kiel_ctdconvert.o:	kiel_ctdconvert.c $(HYDRO_H) 
	$(CC) $(CFLAGS) $(DEBUG) kiel_ctdconvert.c
#
sbe_ctdconvert.o:	sbe_ctdconvert.c $(HYDRO_H) 
	$(CC) $(CFLAGS) $(DEBUG) sbe_ctdconvert.c
#
mmp_convert.o:  mmp_convert.c  $(HYDRO_H)  $(MEMORY_H)
	$(CC) $(CFLAGS) $(DEBUG) mmp_convert.c
#
odf_ctdconvert.o:	odf_ctdconvert.c $(HYDRO_H) $(FILTERS_H)
	$(CC) $(CFLAGS) $(DEBUG) odf_ctdconvert.c
#
pe_transport.o:	pe_transport.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) pe_transport.c
#
soc_ctdconvert.o:	soc_ctdconvert.c $(HYDRO_H)  $(MEMORY_H)
	$(CC) $(CFLAGS) $(DEBUG) soc_ctdconvert.c
yearsmooth.o:	yearsmooth.c
	$(CC) $(CFLAGS) $(DEBUG) yearsmooth.c
#
#  QC PROGRAMS
#
hb_rangechk_ts:	hb_rangechk_ts.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_rangechk_ts.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_rangechk_ts $(LIBS)
#
hb_statchk_ts:	hb_statchk_ts.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_statchk_ts.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_statchk_ts $(LIBS)
#
hb_statfit_ts:	hb_statfit_ts.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_statfit_ts.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_statfit_ts $(LIBS)
#
hb_rangechk_to:	hb_rangechk_to.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_rangechk_to.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_rangechk_to $(LIBS)
#
hb_statchk_to:	hb_statchk_to.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_statchk_to.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_statchk_to $(LIBS)
#
hb_statfit_to:	hb_statfit_to.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS) hb_statfit_to.o hydro_utils.o prop_subs.o phyprops.o eos80.o -o hb_statfit_to $(LIBS)
#
woce_btl2hb_qc:	woce_btl2hb_qc.o hydro_utils.o prop_subs.o phyprops.o eos80.o
	$(LINK) $(LDFLAGS)  woce_btl2hb_qc.o phyprops.o prop_subs.o hydro_utils.o eos80.o $(CDF_L) -o woce_btl2hb_qc $(LIBS)
#
woce_btl_qualchk:	woce_btl_qualchk.o 
	$(LINK) $(LDFLAGS)  woce_btl_qualchk.o $(CDF_L) -o woce_btl_qualchk $(LIBS)
#
woce_btl_qtally:	woce_btl_qtally.o 
	$(LINK) $(LDFLAGS)  woce_btl_qtally.o $(CDF_L) -o woce_btl_qtally $(LIBS)
#
#
#  QC OBJECT DEPENDENCIES
#
hb_rangechk_ts.o:	hb_rangechk_ts.c $(HYDRO_H) 
	$(CC)  $(CFLAGS) $(DEBUG) hb_rangechk_ts.c
#
hb_statchk_ts.o:	hb_statchk_ts.c $(HYDRO_H) $(SIGBINS_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_statchk_ts.c
#
hb_statfit_ts.o:	hb_statfit_ts.c $(HYDRO_H) $(SIGBINS_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_statfit_ts.c
#
hb_rangechk_to.o:	hb_rangechk_to.c $(HYDRO_H) 
	$(CC)  $(CFLAGS) $(DEBUG) hb_rangechk_to.c
#
hb_statchk_to.o:	hb_statchk_to.c $(HYDRO_H) $(SIGBINS_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_statchk_to.c
#
hb_statfit_to.o:	hb_statfit_to.c $(HYDRO_H) $(SIGBINS_H)
	$(CC)  $(CFLAGS) $(DEBUG) hb_statfit_to.c
#
woce_btl2hb_qc.o:	woce_btl2hb_qc.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) woce_btl2hb_qc.c
#
woce_btl_qualchk.o:	woce_btl_qualchk.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) woce_btl_qualchk.c
#
woce_btl_qtally.o:	woce_btl_qtally.c $(HYDRO_H)
	$(CC)  $(CFLAGS) $(DEBUG) woce_btl_qtally.c
#
