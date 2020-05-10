%% Generate the coverage matrix used for the coverage constraint.
%
% Args:
%   O: list of PoIs to monitor, [x, y] in each row
%   N: the struct to record the grid space
%   S_r: the binary sensing range
%
% Return:
%   A: the generated N_o x N_cnt matrix for coverage constraint
%      A * x >= k * 1

function A = cover_matrix(O, N, S_r)

N_o = size(O, 1);       % number of PoIs
N_cnt = size(N, 1);     % number of grid points
A = zeros(N_o, N_cnt);  % initialize A

for i = 1:N_o
    for j = 1:N_cnt
        if norm(O(i) - N(j).position) < S_r
            A(i, j) = 1;
        end
    end
end
end