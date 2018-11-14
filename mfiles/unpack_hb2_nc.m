function [P, T, S, Np, Nt, Ns, lon, lat, stddepth,gridinfo] = unpack_hb2_nc(fname, nan_missing, nan_masked)
%function [P, T, S, Np, Nt, Ns,lon,lat, stddepth,gridinfo] = unpack_hb2_nc(fname, nan_missing, nan_masked);
%OR  function [P, T, S, Np, Nt, Ns,lon,lat, stddepth,gridinfo] = unpack_hb2_nc(fname);
%
% Loads gridded fields from HydroBase 2 netcdf file
%      fname:  name of nc file
%   nan_missing, nan_masked: if  non-zero, sets flagged nodes to NaN
%   if these args are set to zero or not specified, output arrays contain mask and missing flags.
%
% OUTPUT:
%  P, T, S:  3-D matrices of pressure, temperature, salinity 
%                     (dimensions:  nz, ny, nx)
%  Np, Nt, Ns :  number of observations for each gridpoint
%  lon, lat, stddepth : coordinate (column) vectors 
% gridinfo:  struct with fields nx, ny, nz,  xmin, xmax, xinc,ymin,ymax, yinc
%  corresponding to grid bounds not grid nodes (lon,lat vectors contain grid node
%  locations).
 P= NaN; T= NaN; S= NaN;
 Np = NaN; Nt=NaN; Ns=NaN; lon=NaN; lat=NaN; stddepth=NaN; gridinfo.nx = 0;
 
if (nargin ~= 1) && (nargin ~=3)
    error('wrong number of input args');
end
    
if (nargin == 1)
    nan_missing = 0;
    nan_masked = 0;
end

if (nargout ~= 10)
    error('Wrong number of output args');
end

ncid = netcdf.open(fname, 'NC_NOWRITE');

node_offset =  netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'node_offset');

% define latitude vector from north to south (latmax is 1st element) 
mi = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'latmin');
ma =  netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'latmax');
incr =netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'latincr');

if node_offset == 1  % for input pixel grid, set lat vector to reflect node offset
   gridinfo.ymin = mi;
   gridinfo.ymax = ma;
   ma = ma - 0.5 * incr;
   mi = mi + 0.5 * incr;
else  % for input node grid, extend ymin/ymax by 1/2 grid inc to make a pixel grid
  gridinfo.ymin = mi - 0.5 * incr;
  gridinfo.ymax = ma + 0.5 * incr;
end

gridinfo.yinc = incr;
lat = cast( (ma:-incr:mi)', 'double');
gridinfo.ny= length(lat);

% longitude vector (west to east)
mi = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'lonmin');
ma =  netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'lonmax');
incr =netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'lonincr');
if node_offset == 1  % for input pixel grid, adjust lon vector to reflect node offset
   gridinfo.xmin = mi;
   gridinfo.xmax = ma;
   ma = ma - 0.5 * incr;
   mi = mi + 0.5 * incr;
else  % for input node grid, extend min/max by 1/2 grid inc to make a pixel grid
  gridinfo.xmin = mi - 0.5 * incr;
  gridinfo.xmax = ma + 0.5 * incr;
end
gridinfo.xinc = incr;
lon = cast((mi:incr:ma)', 'double');
gridinfo.nx= length(lon);
clear mi ma incr

gridinfo.node_offset = 1;  % standardize this for all hb grids

% get depth dimension 
dimid_de = netcdf.inqDimID(ncid,'de');
[~, nlevs] = netcdf.inqDim(ncid,dimid_de);
gridinfo.nz = nlevs;
varid = netcdf.inqVarID(ncid,'de');
X=  netcdf.getVar(ncid,varid);
stddepth=cast(X,'double');

% 3D grids are stored in the order :  nz, nx, ny
% Permute the dimensions to nz, ny, nx

varid = netcdf.inqVarID(ncid,'pr');
X=  netcdf.getVar(ncid,varid);
P=cast(permute(X,[1 3 2]), 'double');
if (nan_masked)
   P(P>10000)=NaN;
end
if (nan_missing)
    P(P<-10000)=NaN;
end

varid = netcdf.inqVarID(ncid,'pr_cnt');
N= netcdf.getVar(ncid,varid);
Np = permute(N,[1 3 2]);

varid = netcdf.inqVarID(ncid,'te');
X= netcdf.getVar(ncid,varid);
T=cast(permute(X,[1 3 2]), 'double');
if (nan_masked)
   T(T>10000)=NaN;
end
if (nan_missing)
    T(T<-10000)=NaN;
end

varid = netcdf.inqVarID(ncid,'te_cnt');
N= netcdf.getVar(ncid,varid);
Nt = permute(N,[1 3 2]);

varid = netcdf.inqVarID(ncid,'sa');
X= netcdf.getVar(ncid,varid);
S=cast(permute(X,[1 3 2]), 'double');
if (nan_masked)
   S(S>10000)=NaN;
end
if (nan_missing)
    S(S < -10000)=NaN;
end
varid = netcdf.inqVarID(ncid,'sa_cnt');
N= netcdf.getVar(ncid,varid);
Ns = permute(N,[1 3 2]);



