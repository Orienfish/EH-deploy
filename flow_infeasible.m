%% Get the equality constraints for setting the infeasible flow 
%  to zero.
% Args:
%   N: the struct of grid locations
%   v_cnt: total number of variables in the optimization problem
%   xform: the format of the variables
%   C_r: allowed communication range in m
%
% Return:
%   Aeq, beq: the equality constraint fed to the solver to filter out
%             infeasible flows. Each constraint set one flow to zero.

function [Aeq, beq] = flow_infeasible(N, v_cnt, xform, C_r)
% get the number of grid points
N_cnt = size(N, 1);

% iteratively construct Aeq line by line
% each line corresponds to one infeasible flow
Aeq = [];
for i = 1:N_cnt
    for j = i+1:N_cnt
        if norm(N(i).position - N(j).position) > C_r
            % set the flow from node i to j to zero
            line = zeros(1, v_cnt);
            line(xform.fij_base + i * N_cnt + j) = 1;
            Aeq = [Aeq; line];
            % set the flow from node j to i to zero
            line = zeros(1, v_cnt);
            line(xform.fij_base + j * N_cnt + i) = 1;
            Aeq = [Aeq; line];
            fprintf('%d, %d\n', i, j);
        end
    end
end

% set beq
beq = zeros(size(Aeq, 1), 1);
end

