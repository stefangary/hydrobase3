function [profile, err] = hb_get_profile(fid)
% Retrieves the next profile from an already opened HydroBase 
% text file.  Returns a structure with the following form:
%  profile.country
%         .ship, 
%         .cruise 
%         .station
%         .year  
%         .month  
%         .day 
%         .orig 
%         .instr
%         .lat  
%         .lon 
%         .bdpth  
%         .nobs 
%         .nprops 
%         .ms10  
%         .ms1 
%         .prop_id  ( names of variable fields in profile.data)
%         .qual 
%         .data  (structure containing nprops separate fields )
%              .pr
%              .t90  (etc....)
% Returns err= 0 for successful read, -1 for EOF.
% If an error occurs, it displays a message and returns err = 1;
%____________________________________________________________
%%  check args and file status
profile = [];
err = 0;

if nargin ~= 1
   error('hb_get_profile():  Must supply already opened file id');
end

if feof(fid)
   err = -1;
   return
end
% 
[profile, err] = hb_read_hdr(fid);
if err
   return;
end
%% Read all data levels into a  cell array and assign columns of cell to
%% fields of struct profile.data.(prop_mne).  Convert HB_MISSING values to NaN.

field = '%f ';
fmt = [];
for ii = 1:profile.nprops
   fmt = [fmt field];
end

C = textscan(fid,fmt, profile.nobs);

for ii=1:profile.nprops
    mne = strtrim(profile.prop_id(ii,:));
    profile.data.(mne) = C{1,ii};
    missing = find(profile.data.(mne) < -8.9999);
    if ~isempty(missing)
        profile.data.(mne)(missing) = NaN;
    end
 
    
end

delim = fgetl(fid);  %move past lf char
delim = fgetl(fid);

if ~ (findstr('**', delim)) 
    err = 1;
    warning('Error reading station delimiter in hb_get_profile()');
end

% test for end of file
if feof(fid)
   err = -1;
end

% End of function
return

