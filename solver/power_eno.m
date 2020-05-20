%% Get the power inequality constraint to satisfy ENO
% Args:
%   eta: the sampling frequency
%   B: bandwidth B/s
%   v_cnt: the total number of variables
%   xform: the format of the variables
%   N: the struct of grid locations
%   dist: the distance matrix between all grid locations and the sink
%   C_r: the communication radius
%   rel: the struct of reliability options and targets
%
% Return:
%   A, b: a N_cnt * v_cnt matrix and a N_cnt * 1 vector for power ENO
%         constraint A * x <= b

function [A, b] = power_eno(eta, B, v_cnt, xform, N, dist, C_r, rel)
% initialization
N_cnt = xform.x_cnt;         % number of grid locations


% constant parameters setting
P0 = 0.01;                  % sleep power (W)
Es = 0.2*0.2;               % sensing power (W*s)
Prx = 0.1;                  % reception power (W)

% get the transformation matrix to extract s
T_s = [zeros(xform.s_cnt, xform.x_cnt), eye(xform.s_cnt), ...
    zeros(xform.s_cnt, v_cnt-xform.s_end)];
% N_cnt * v_cnt matrix for si*eta*Es
A = eta * Es * T_s;

% set the coefficients for the flow part
for i = 1:N_cnt
    for j = 1:N_cnt
        if i == j % avoid the same node
            continue;
        end
        if dist(i, j) <= C_r
            % outgoing flow from node i to node j, fij
            fij_idx = xform.fij_base + (i-1)*N_cnt + j;
            A(i, fij_idx) = (1/B) * getPtx(dist(i, j));
            % incoming flow from node j to node i, fji
            fji_idx = xform.fij_base + (j-1)*N_cnt + i;
            A(i, fji_idx) = (1/B) * Prx;
        end
    end
    if dist(i, N_cnt+1) <= C_r
        % outgoing flow from node i to the sink, fiB
        fiB_idx = xform.fiB_base + i;
        A(i, fiB_idx) = (1/B) * getPtx(dist(i, N_cnt+1));
    end
end
%fprintf('power eno:\n');
%disp(A);
b = vertcat(N(:).Ri) - repmat(P0, N_cnt, 1);
if rel.SoH == true
    P_soh = Psoh_bound(rel.SoHref, rel.T, vertcat(N(:).Ti));
    b = [b, P_soh - repmat(P0, N_cnt, 1)];
    b = max(b, 0);
end
if rel.MTTF == true
    P_mttf = Pmttf_bound(rel.MTTFref, vertcat(N(:).Ti));
    b = [b, P_mttf - repmat(P0, N_cnt, 1)];
    b = max(b, 0);
end
disp(b);
b = min(b, [], 2); % get the column vector of min of each row
end

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


