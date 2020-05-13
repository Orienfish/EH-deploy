%% Get the flow bound inequality constraint. 
%  Ensure that the total amount of flow is zero at undeployed locations
% Args:
%   v_cnt: the total number of variables
%   xform: the format of the variables
%   dist: the distance matrix of all grid locations and the sink
%   C_r: the feasible communicaition radius
%   maxf: the maximum flow amount
%
% Return:
%   A, b: the generated N_cnt * v_cnt matrix and N_cnt * 1 vector for
%         flow conservation constraint.

function [A, b] = flow_bound(v_cnt, xform, dist, C_r, maxf)
% initialization
N_cnt = xform.x_cnt;           % number of grid nodes
A = zeros(N_cnt, v_cnt);

% set the deployed/undeployed part
for i = 1:N_cnt
    A(i, i) = -maxf;
end

% set the coefficients for the flow part
for i = 1:N_cnt
    for j = 1:N_cnt
        if i == j % avoid the same node
            continue;
        end
        if dist(i, j) <= C_r
            % outgoing flow from node i to node j, fij
            fij_idx = xform.fij_base + (i-1)*N_cnt + j;
            A(i, fij_idx) = 1;
        end
    end
    if dist(i, N_cnt+1) <= C_r
        % outgoing flow from node i to the sink, fiB
        fiB_idx = xform.fiB_base + i;
        A(i, fiB_idx) = 1;
    end
end
%fprintf('flow bound:\n');
%disp(A);
b = zeros(N_cnt, 1);
end

