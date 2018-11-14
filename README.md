# hydrobase3
A data processing toolkit for oceanographers

This repository started with the Hydrobase 3 files
of the HB3_pkg.2013may22.tar.gz distribution on the
WHOI website.  Hydrobase can be distributed via
GitHub provided that it remains open source.

The only files from that distribution that
I did not include are ./lib/etopo1_ice_gline.grd
and ./lib/global_oi_parms.nc because they are
each more than 100MB.

The former can be reconstructed by
1) downloading ETOPO1_Ice_g_gmt4.grd.gz from
https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/
2) and then running the GMT command:
grdreformat ETOPO1_Ice_g_gmt4.grd=ns/1/0/-32768 /data1/HB3/lib/etopo1ice__gline.grd=ns/1/0/-32768

