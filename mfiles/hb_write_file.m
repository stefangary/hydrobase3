function err = hb_write_file(profiles, fname, mode)
% function err = hb_write_file(profiles, fname, mode)
%
% Creates file for writing in specified mode, writes profiles
% and closes the file.
%  profiles is struct array of profiles to be written
%  fname is full pathname of file
%  mode is 0 (zero) or 'w' for Overwrite
%          1 NoClobber
%          2 or 'a' for Append
%         
% Returns err > 0 for an error.
err = 1;
if nargin < 2
    help hb_write_file
    error('Fatal error')
end

if nargin == 2
    mode = 1;
end

fidout = hb_create_file(fname, mode);
[m,nsta] = size(profiles);
for ii=1:nsta
    err = hb_write_profile(fidout,profiles(ii));
end
disp(['finished writing ' fname]);

fclose(fidout);



