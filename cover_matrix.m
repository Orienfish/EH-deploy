%% Generate the coverage matrix used for the coverage constraint.
%
% Args:
%   O: list of PoIs to monitor, [x, y] in each row
%   N: the struct to record the grid space
%   S_r: the binary sensing range
%   v_cnt: the total number of variables
%
% Return:
%   A: the generated N_o * v_cnt matrix for coverage constraint.
%      Note that only the first N_cnt columns matter in the coverage constraint.
%      A * x >= k * 1

function A = cover_matrix(O, N, S_r, v_cnt)

N_o = size(O, 1);       % number of PoIs
N_cnt = size(N, 1);     % number of grid points
A = zeros(N_o, v_cnt);  % initialize A

for i = 1:N_o
    for j = 1:N_cnt
        if norm(O(i, :) - N(j).position) <= S_r
            %fprintf('O %d: (%f %f)\n', i, O(i, 1), O(i, 2));
            %fprintf('N %d: (%f %f)\n', j, N(j).position(1), N(j).position(2));
            %fprintf('dist: %f\n', norm(O(i, :) - N(j).position));
            A(i, j) = 1;
        end
    end
end
end