function [index] = hb_get_prop_indx(str)
%%

if isempty(who('HB_INITIALIZED', 'global'))
    hb_initialize;
end
if isempty(who('HB_INITIALIZED'))
    global HB_INITIALIZED HB_PROP_INFO
end
%%

if (nargin ~= 1) | (~isstr(str))
    error('hb_find_prop_indx(): Must supply property name as a string')
end
%%

index=0;
[m,n] = size (HB_PROP_INFO);
for ii=1:m
    if strcmp(str, HB_PROP_INFO(ii,1));
       index = ii;
       return
    end
end

return
% End of function
