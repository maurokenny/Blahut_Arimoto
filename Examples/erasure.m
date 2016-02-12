clear p;
addpath ../
% For example, the transition matrix for the erasure channel is
% can be calculated as
e = 0.5;
p = [1-e e 0; 0 e 1-e]; % conditional prob. for erasure channel
% The capacity can be calculated by BlahutArimoto, and is equal to 1-e
BlahutArimotoConstraint(p,[0 1],1,2)

