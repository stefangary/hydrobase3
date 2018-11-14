function fid = hb_create_file(fname, mode)
%%%  Opens a file for writing in specified mode.
%%%  mode = 
%%%         0 = overwrite if file exists
%%%         1 = no_clobber (do not overwrite if file exists)
%%%         2 = append or create
%%%        'w' = overwrite or create
%%%        'a' = append or create
%%%  An error will occur if mode is 1 (noclobber) and the file exists.
%%%  Returns a valid file id, or -1 for an error.
%%%
%%% Usage:   fid = open_hydro_file(fname, mode)
%%  check args
fid = -1;
if nargin ~= 2
    help hb_create_file
    error('Error in call to hb_create_file()');
end

%%  create file
switch (mode)
    case 0      
        [fid, msg] = fopen(fname, 'w');
         if fid == -1
             error(msg)
         end
    case 1
        % Test for existence of file
        [fid, msg] = fopen(fname, 'r');
        if fid >= 0
            fclose(fid);
            error([ 'FATAL ERROR: ', fname, ' already exists.']);
        end
        [fid, msg] = fopen(fname, 'w');
         if fid == -1
             error(msg)
         end
    case 2
        
        [fid, msg] = fopen(fname, 'a');
         if fid == -1
             error(msg)
         end
    case {'w','a','w+','a+'}
        [fid, msg] = fopen(fname, mode);
         if fid == -1
             error(msg)
         end
    otherwise
        disp('Invalid file mode passed to hb_create_file()')
end

