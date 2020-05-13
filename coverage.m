%% Generate the coverage inequality constraint.
%
% Args:
%   O: list of PoIs to monitor, [x, y] in each row
%   N: the struct to record the grid space
%   S_r: the binary sensing range
%   v_cnt: the total number of variables
%   K: predefined coverage level
%
% Return:
%   A, b: the generated N_o * v_cnt matrix and N_o * 1 vector for coverage 
%         constraint.
%         Note that only the first N_cnt columns matter in the coverage 
%         constraint.
%         A * x >= k * 1

function [A, b] = coverage(O, N, S_r, v_cnt, K)

N_o = size(O, 1);       % number of PoIs
N_cnt = size(N, 1);     % number of grid points
A = zeros(N_o, v_cnt);  % initialize A

for i = 1:N_o
    for j = 1:N_cnt
        if norm(O(i, :) - N(j).position) <= S_r
            %fprintf('O %d: (%f %f)\n', i, O(i, 1), O(i, 2));
            %fprintf('N %d: (%f %f)\n', j, N(j).position(1), N(j).position(2));
            %fprintf('dist: %f\n', norm(O(i, :) - N(j).position));
            s_idx = N_cnt + j;
            A(i, s_idx) = 1;
        end
    end
end
disp(A);
A = -A;
b = -repmat(K, N_o, 1);
end