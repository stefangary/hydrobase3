function [p, err] = hb_read_hdr(fid)
% Reads header info from an already opened HydroBase station file
% and store it in fields of struct p.  
% Returns err = 0 for a successful read, -1 for EOF.
% If an error occurs, writes a message and returns 1.
%%  Check input
p = [];
err = 0;

if nargin ~= 1
   error('hb_read_hdr():  Must supply already opened file id');
end

if feof(fid)
   err = -1;
   return
end
%% Read first header line

buf = fgetl(fid);
[m,n] = size( buf );
if n~= 75
   warning(['Not a valid HydroBase header record:' char(buf{1})] )
   err = 1;
   return
end

H = textscan(buf,'%s %s %d %d %d %d %d %c%c %f %f %d %s %d %d %d %d',1);
[m,n] = size(H);

p.country = char(H{1});
p.ship = char(H{2});
p.cruise = H{3};
p.station = H{4};
p.year = H{5};
p.month  = H{6};
p.day  = H{7};

if p.year < 1000
    p.year = p.year + 1900;
end

if ~isempty(H{17})   %  no missing fields, parse the remainder
    p.orig  = char(H{8});
    p.instr  = char(H{9});
    p.lat  = H{10};
    p.lon  = H{11};
    p.bdpth  = H{12};
    p.qual  = char(H{13});
    p.nobs  = H{14};
    p.nprops  = H{15};
    p.ms10  = H{16};
    p.ms1  = H{17};
else
    line = char(buf{1});
    p.ms10 = sscanf(line(69:72),'%d');
    
    if p.ms10 > 1000  % HydroBase2 header with 1 or more blank char fields
        p.orig = sscanf(line(29),'%c');
        p.instr = sscanf(line(30),'%c');
        a = sscanf(line(32:53),'%f %f %d');
        p.lat  = a(1);
        p.lon  = a(2);
        p.bdpth  = a(3);
        a = sscanf(line(60:75),'%d %d %d %d');
        p.nobs = a(1);
        p.nprops = a(2);
        p.ms10 = a(3);
        p.ms1 = a(4);
        p.qual = sscanf(line(55:58),'%4s');
        
    else  % Old HydroBase header
        a = sscanf(line(31:75),'%f %f %d %d %d %d %d');
        p.lat = a(1);
        p.lon = a(2);
        p.bdpth = a(3);
        p.nobs = a(4);
        p.nprops = a(5);
        p.ms10 = a(6);
        p.ms1 = a(7);
        p.orig = '0';
        p.instr = 'u';
        p.qual = '0000';
     end
end

clear H, buf;
%%  Read line 2: the list of properties

buf = fgetl(fid);
H = textscan(buf, '%s');

[nprops, ncols] = size(H{1});
if (nprops ~= p.nprops)
    warning('hb_read_hdr(): Mismatch between nprops field in hdr line1 and # of props listed in hdr line2');
    err = 1;
    return
end
p.prop_id = char(H{1});


% End of function
return

    