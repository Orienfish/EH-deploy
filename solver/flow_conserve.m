%% Get the flow conservation equality constraint
% Args:
%   G: the generated bytes of each sample
%   eta: the sampling frequency
%   v_cnt: the total number of variables
%   xform: the format of the variables
%   dist: the distance matrix of all grid locations and the sink
%   C_r: the feasible communicaition radius
%
% Return:
%   A, b: the generated N_cnt * v_cnt matrix and N_cnt * 1 vector for
%         flow conservation constraint A * x >= b
function [A, b] = flow_conserve(G, eta, v_cnt, xform, dist, C_r)
% initialization
N_cnt = xform.x_cnt;          % number of grid nodes

% get the transformation matrix to extract s
T_s = [zeros(xform.s_cnt, xform.x_cnt), eye(xform.s_cnt), ...
    zeros(xform.s_cnt, v_cnt-xform.s_end)];
% N_cnt * v_cnt matrix for si*eta*G
A = G * eta * T_s;

% set the coefficients for the flow part
for i = 1:N_cnt
    for j = 1:N_cnt
        if dist(i, j) <= C_r
            % outgoing flow from node i to node j, fij
            fij_idx = xform.fij_base + (i-1)*N_cnt + j;
            A(i, fij_idx) = -1;
            % incoming flow from node j to node i, fji
            fji_idx = xform.fij_base + (j-1)*N_cnt + i;
            A(i, fji_idx) = A(i, fji_idx) + 1; % so that fii = 0
        end
    end
    if dist(i, N_cnt+1) <= C_r
        % outgoing flow from node i to the sink, fiB
        fiB_idx = xform.fiB_base + i;
        A(i, fiB_idx) = -1;
    end
end
%fprintf('flow conservation:\n');
%disp(A);
b = zeros(N_cnt, 1);
end

