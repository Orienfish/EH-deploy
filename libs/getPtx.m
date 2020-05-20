%% Get the transmission power from distance
% Args:
%   dist_m: the distance in m
%
% Return:
%   Ptx: the transmission power in W

function Ptx = getPtx(dist_m)
Pto = 0.22;                 % transmission power baseline
alpha = 3.5;                % distance coefficient
beta = 1e-7;                % distance coefficient

Ptx = Pto + beta * dist_m ^ alpha;
end