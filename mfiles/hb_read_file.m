function [profiles,nsta] = hb_read_file(dir,root,ext )
% function [profiles,nsta] = hb_read_file(dir,root,ext )
% Reads a complete HydroBase ascii station file into a struct array of profiles.
% Filename is specified EITHER by 3 strings for directory, root, and extent (which 
% are concatenated to form a (platform dependent) filename OR by a single string 
% specifying the entire filename.
%
% Returns profiles as a struct array and the number of stations read in.
% USAGE [p,n] = hb_read_file('d1/HB2/Data', '5603', '.btl');
%    OR [p,n] = hb_read_file('d1/HB2/Data/5603.btl');


profiles = [];
nsta = 0;
if nargin ~= 1 & nargin ~= 3
    help hb_read_file
    error('Accepts filename as either 3 args or 1')
end

if nargin == 1
    fidin = hb_open_file(dir);
else
    fidin = hb_open_file(dir, root, ext);
end

if fidin < 0
    return
end
err = 0;
[profiles,err] = hb_get_profile(fidin);
nsta = 1;
while (err ~= -1)
    if err > 0
        disp('FATAL ERROR in hb_read_file');
        error(['Stopped after ' num2str(nsta) ' stations read in'])
    end
    nsta = nsta + 1;
    [profiles(nsta),err] = hb_get_profile(fidin);
end

fclose(fidin);



