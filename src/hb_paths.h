/* hb_paths.h
................................................................................
                          *******  HydroBase 3 *******
                             Ruth Curry
                             Woods Hole Oceanographic Institution
			     updated May 2010
			 
...................................................

2 versions of ETOPO1 topography files are supplied with HydroBase3:
 -- ice surface versions, grid-registered -- 
 full resolution (1 min) and coarse resolution versions (0.1-degree).   
  
GAMMA_NC_PATH is a data file used for computing neutral density (gamma_n)

LISTPATH is the directory for HydroBase files which define isopycnal ranges for each 10-deg square

OI_PARMS_FILE contains variogram and gradient parameters for use in optimal interpolation.
*/

#define LIBPATH "/usr/local/hb3/lib/"
#define LISTPATH "/usr/local/hb3/lists/"

#define BATHPATH "/usr/local/hb3/lib/etopo1_ice_gline.grd"
#define BATHPATH_C "/usr/local/hb3/lib/topo_tenthdeg_ice_gline.grd"

#define GAMMA_NC_PATH "/usr/local/hb3/lib/gamma.nc"
#define SHIPCODE_PATH "/usr/local/hb3/lib/wod_shipcodes.txt"

#define OI_PARMS_FILE "/usr/local/hb3/lib/global_oi_parms.nc"
