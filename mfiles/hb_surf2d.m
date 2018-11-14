function data = hb_surf2d(infiles, surfprop, surfval, otherargs)
% function data = hb_surf2d(infiles, prop, val, otherargs)
%   Usage:  hb_surf2d filename_root(s)  [-D<dirname>] [-E<file_extent>] 
%              -P<list_of_properties>  [-Y] [-M[d]]  [-L]   [-I] 
%             [-B/west/east/south/north] [-S<surface_def_file>] [-W<window>[/<w_incr>]]
%           
%   infiles :  string containing ms10list,  -D and  -E info
%  surfprop : a string containing prop mnemonic (and pref if appropriate)
%  surfval : scalar 
%  otherargs :  string containing -P<proplist> and other options 

if nargin < 4
    system('hb_surf2d -h');
    error('USAGE:  hb_surf2d(infiles, surfprop, surfval, otherargs) ');
end

fid = fopen('surf.lis','wt');
fprintf(fid,'%s  %f %s \nend\n', surfprop, surfval, 'hb_surf2d.out');
fclose(fid);

string = ['hb_surf2d ', infiles, '   -Ssurf.lis -L ', otherargs];
system(string);

data = load('hb_surf2d.out');