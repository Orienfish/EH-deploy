%% Get the power inequality constraint to satisfy ENO
% Args:
%   v_cnt: the total number of variables
%   xform: the format of the variables
%   N: the struct of grid locations
%   dist: the distance matrix between all grid locations and the sink
%   params: necessary basic parameters
%   rel: the struct of reliability options and targets
%
% Return:
%   A, b: a N_cnt * v_cnt matrix and a N_cnt * 1 vector for power ENO
%         constraint A * x <= b

function [A, b] = power_eno(v_cnt, xform, N, dist, params, rel)
% initialization
N_cnt = xform.x_cnt;         % number of grid locations

% get the transformation matrix to extract s
T_s = [zeros(xform.s_cnt, xform.x_cnt), eye(xform.s_cnt), ...
    zeros(xform.s_cnt, v_cnt-xform.s_end)];
% N_cnt * v_cnt matrix for si*eta*Es
A = params.eta * params.Es * T_s;

% set the coefficients for the flow part
for i = 1:N_cnt
    for j = 1:N_cnt
        if i == j % avoid the same node
            continue;
        end
        if dist(i, j) <= params.C_r
            % outgoing flow from node i to node j, fij
            fij_idx = xform.fij_base + (i-1)*N_cnt + j;
            A(i, fij_idx) = (1/params.B) * getPtx(dist(i, j));
            % incoming flow from node j to node i, fji
            fji_idx = xform.fij_base + (j-1)*N_cnt + i;
            A(i, fji_idx) = (1/params.B) * params.Prx;
        end
    end
    if dist(i, N_cnt+1) <= params.C_r
        % outgoing flow from node i to the sink, fiB
        fiB_idx = xform.fiB_base + i;
        A(i, fiB_idx) = (1/params.B) * getPtx(dist(i, N_cnt+1));
    end
end
%fprintf('power eno:\n');
%disp(A);
if rel.SoH && rel.MTTF        % conside reliability constraints
    b = vertcat(N(:).Pi) - repmat(params.P0, N_cnt, 1);
else                          % only energy neutral operation constraint
    b = vertcat(N(:).Ri) - repmat(params.P0, N_cnt, 1);
end
%disp(b);
end


