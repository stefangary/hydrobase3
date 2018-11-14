function status = hb_write_hdr(fid, p)
%%% Writes header info to an already opened HydroBase3 station
%%% file.  Returns 0 for a successful write or a value > 0 if an error
%%% occurs.
%%%
%%% Usage:  status = hb_write_hdr(fid, p)

status = 0;
if nargin < 2
        help hb_write_hdr
        status = 1;
    error('Error in call to hb_write_hdr()');
end

if ~(isstruct(p) & isstruct(p.data))
    status = 1;
    warning('Invalid profile structure passed to hb_write_hdr()');
    return
end
%%  check that year is 4 digits (not 2)
if p.year < 1000
    p.year = p.year + 1900;
end
%% check quality bytes
if ~ischar(p.qual)
    p.qual = ['0', '0', '0','0'];
else
   for ii=1:4
       if (p.qual(ii) < '0' | p.qual(ii) > '9' | p.qual == NaN )
        p.qual(ii) = '0';
       end 
   end
end
%% check origin and instrument fields
if ~ischar( p.orig)
    p.orig = '0';
end

if ~ischar(p.instr)
    p.instr = 'u';
end
%% compose first header line and write to file

buf = sprintf('%.2s %.2s %5d %4d %4d %2d %2d %c%c %7.3f %8.3f %5d %4s %4d %3d %4d %2d\n',...
    p.country, p.ship, p.cruise, p.station, p.year, p.month, p.day, p.orig, p.instr, p.lat,...
    p.lon, p.bdpth, p.qual, p.nobs, p.nprops, p.ms10, p.ms1);

if length(buf) ~= 76
    status = 1;
    error('Fatal Error writing first header line of profile (incorrect # of bytes)');
end

%% compose line of property descriptors and write to file...
buf = [buf, p.prop_id(1,:)];

for ii=2:p.nprops
buf = [buf, ' ', p.prop_id(ii,:)];
end

str = sprintf('\n');
buf = [buf,str];
n = fwrite(fid,buf,'uchar');

return
% End of function




