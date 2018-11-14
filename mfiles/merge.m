function out = merge(input1, input2)
%function out = merge(input1, input2)
% combines two matrices size by side and pads the one with smaller length
% with NaNs at the bottom 

%P Robbins 2/95
out=[]; l1 = []; l2 = []; w1 = []; w2 = [];

[l1,w1] = size(input1);
[l2,w2] = size(input2);


  if l1 == l2
    out = [input1 input2];
  elseif l1 > l2
    out = [input1 [input2; nan*ones(l1-l2,w2)]];
  else
    out = [[ input1 ;nan*ones(l2-l1,w1)] input2];
  end
