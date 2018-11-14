function [ms10, ms1] = hb_ms10(lat, lon)
% Returns the 10-deg and 1-deg WMO square designations for a particular
% latitude, longitude position.
%
%  USAGE:  [ms10, ms1] = hb_ms10(lat, lon)
%% check input args
ms10 = -1;
ms1 = -1;
if nargin < 2
    help hb_ms10
    error('Error in call to hb_ms10()')
end
if (lat > 90) | (lat < -90)
    error('Fatal error in hb_ms10: Latitude not in range -90 to +90');
end
if (lon < -180) | (lon > 360)
    error('Fatal error in hb_ms10: Longitude not in range -180 to +360');
end

%% determine earth quadrant
quadrant = [7000, 1000; 5000, 3000];
k = 1;
if lat < 0
    k = 2;
end
kk = 2;
if lon < 0
    kk = 1;
end
if lon >= 180.
    lon = lon - 360;
    kk = 1;
end
%% this smidgeon handles borderline cases correctly in all hemispheres
ilat = floor(lat + .00001);
ilon = floor(lon + .00001);

if ilat < 0
    ilat = -ilat;
end
if ilon < 0
    ilon = -ilon;
end
if ilat >= 90  %special case at the poles
    ilat = 89;
end
ms1 = mod(ilat,10) * 10 + mod(ilon, 10);
ms10 = quadrant(k,kk) + fix(ilat / 10) * 100 + fix(ilon / 10);
return
%% end of function

