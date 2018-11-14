function [P, T, S, Np, Nt, Ns, Ep, Et, Es, lon, lat, stddepth, gridinfo] = unpack_hb3_nc(fname, nan_missing, nan_masked)
%function [P, T, S, Np, Nt, Ns, Ep, Et, Es,lon,lat, stddepth, gridinfo] = unpack_hb3_nc(fname, nan_missing, nan_masked);
%  OR function [P, T, S, Np, Nt, Ns, Ep, Et, Es,lon,lat, stddepth, gridinfo] = unpack_hb3_nc(fname);
%
% Loads gridded fields from HydroBase 3 nc file
%      fname:  name of nc file
%   
%   nan_missing, nan_masked: if  non-zero, sets flagged nodes to NaN
%   if these args are set to zero or not specified, output arrays contain mask and missing flags.
%
% OUTPUT:
%  P, T, S:  3-D matrices of pressure, temperature, salinity (dimensions: nz, ny,nx)
%  Np, Nt, Ns :  number of observations for each gridpoint
%  Ep, Et, Ns :  Error for each gridpoint
%  lon, lat, stddepth : coordinate (column) vectors  
%  gridinfo : struct with fields nx, ny, nz,xmin, xmax, xinc,ymin,ymax,yinc
%  corresponding to grid bounds not grid nodes (lon,lat vectors contain grid node
%  locations).

P= NaN; T= NaN; S= NaN;
 Np = NaN; Nt=NaN; Ns=NaN; lon=NaN; lat=NaN; stddepth=NaN; gridinfo.nx = 0;
if (nargin ~= 1 && nargin ~=3)
    error('Wrong number of input args');
end
    
if (nargin == 1)
    nan_missing = 0;
    nan_masked = 0;
end

if (nargout ~= 13)
    error('Wrong number of output args');
end

ncid = netcdf.open(fname, 'NC_NOWRITE');
node_offset =  netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'node_offset');

% define latitude vector from north to south 
varid = netcdf.inqVarID(ncid,'latitude');
X= netcdf.getVar(ncid,varid);
lat = cast(X, 'double');
gridinfo.yinc =netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'latincr');
gridinfo.ymin = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'latmin');
gridinfo.ymax =  netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'latmax');
if node_offset == 0   % for node grids extend bounds by 1/2 gridincrement
    gridinfo.ymin = gridinfo.ymin - 0.5 * gridinfo.yinc;
    gridinfo.ymax = gridinfo.ymax + 0.5 * gridinfo.yinc;
end
gridinfo.ny = length(lat);
gridinfo.node_offset = 1;  % standardize to pixel grid for all hb grids

% longitude vector
varid = netcdf.inqVarID(ncid,'longitude');
X= netcdf.getVar(ncid,varid);
lon = cast(X, 'double');
gridinfo.xinc =netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'lonincr');
gridinfo.nx = length(lon);
gridinfo.xmin = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'lonmin');
gridinfo.xmax =  netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'lonmax');
if node_offset == 0   % for input node grids extend bounds by 1/2 gridincrement
    gridinfo.xmin = gridinfo.xmin - 0.5 * gridinfo.xinc;
    gridinfo.xmax = gridinfo.xmax + 0.5 * gridinfo.xinc;
end
clear X

% get depth dimension 
dimid_de = netcdf.inqDimID(ncid,'de');
[~, nlevs] = netcdf.inqDim(ncid,dimid_de);
gridinfo.nz = nlevs;
varid = netcdf.inqVarID(ncid,'de');
X= netcdf.getVar(ncid,varid);
stddepth = cast(X,'double');

%  we now know the dimensions of each 3D grid which are stored in 
% the order :  nlevs x nlons x nlats
varid = netcdf.inqVarID(ncid,'pr');
X= netcdf.getVar(ncid,varid);
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
varid = netcdf.inqVarID(ncid,'pr_err');
X= netcdf.getVar(ncid,varid);
Ep=cast(permute(X,[1 3 2]), 'double');
if (nan_masked)
   Ep(Ep>10000)=NaN;
end
if (nan_missing)
    Ep(Ep < -10000)=NaN;
end

tname = 't90';
varid = netcdf.inqVarID(ncid,tname);
if (varid < 0)
tname = 'te';
    varid = netcdf.inqVarID(ncid,tname);
end
X= netcdf.getVar(ncid,varid);
T=cast(permute(X,[1 3 2]), 'double');
if (nan_masked)
   T(T>10000)=NaN;
end
if (nan_missing)
    T(T<-10000)=NaN;
end

varid = netcdf.inqVarID(ncid,[tname '_cnt']);
N= netcdf.getVar(ncid,varid);
Nt = permute(N,[1 3 2]);
varid = netcdf.inqVarID(ncid,[tname '_err']);
X= netcdf.getVar(ncid,varid);
Et=cast(permute(X,[1 3 2]), 'double');
if (nan_masked)
   Et(Et>10000)=NaN;
end
if (nan_missing)
    Et(Et < -10000)=NaN;
end

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
varid = netcdf.inqVarID(ncid,'sa_err');
X= netcdf.getVar(ncid,varid);
Es=cast(permute(X,[1 3 2]), 'double');
if (nan_masked)
   Es(Es > 10000)=NaN;
end
if (nan_missing)
    Es(Es < -10000)=NaN;
end


