function fid = hb_open_file(dir, root, ext)
%%%  Opens an existing HydroBase file specified EITHER by 3 character
% strings for directory, root, and extent which are concatenated to form a
% (platform dependent) filename OR by a single string specifying the entire filename
% 
%%%  Returns fid or -1 if not successful.
%%%
%%% Usage:   fid = hb_open_file(dir, root, extent)
%%%   OR:    fid = hb_open_file(filename)
%%
fid = -1;

if nargin ~= 1 & nargin ~= 3
    help hb_open_file
    error('Wrong number of arguments')
end

fname = dir;
if nargin == 3
    fname = fullfile(dir, [root ext]);
end
[fid, msg] = fopen(fname, 'r');
if fid == -1
    disp([msg ' ' fname])
else
    disp(['Opened ' fname ' ...'])
end
