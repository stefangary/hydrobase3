function status = hb_write_profile(fid, p)
%%% Writes one complete profile to an already opened HydroBase3 station
%%% file.  Returns 0 for a successful write or a value > 0 if an error
%%% occurs.
%%%
%%% Usage:  status = hb_write_profile(fid, p)
%%  check for HydroBase global vars

if isempty(who('HB_INITIALIZED', 'global'))
    hb_initialize;
end
if isempty(who('HB_INITIALIZED'))
    global HB_INITIALIZED HB_PROP_INFO HB_MISSING 
end
%% check args
status = 0;
if nargin < 2
        help hb_write_profile
        status = 1;
    error('Error in call to hb_write_profile()');
end

if ~(isstruct(p) & isstruct(p.data))
    status = 1;
    error('Invalid profile structure passed to hb_write_profile()')
end
%%  write header info to file

status = hb_write_hdr(fid, p);
if status > 0
    return
end

%% check for NaNs and replace with -9.0

clear fmt;
clear indx;
for ii=1:p.nprops
    mne = strtrim(p.prop_id(ii,:));
    indx(ii) = hb_get_prop_indx(mne);
    eval(['indnan = find(isnan(p.data.' mne '));'])
    if ~isempty(indnan)
       eval(['p.data.' mne '(indnan) = HB_MISSING;'])
    end
end


%% arrange data row-by-row in formatted buffer
  
  propformat = [' '];
  propvals = repmat(nan,p.nobs,p.nprops);
  for icol=1:p.nprops
    mne = p.prop_id(icol,:);
    propformat = [propformat char(HB_PROP_INFO{indx(icol),4}) ' '];
    propvals(:,icol) = eval(['p.data.' mne ';']);
  end
  propformat = [propformat '\n'];

  fprintf(fid,propformat,propvals');
  fprintf(fid,'**\n');


return

