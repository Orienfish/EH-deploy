%% Get the power at all grid locations
% Args:
%   sol: the struct of the solution to be evaluated
%   N: the struct of grid locations
%   dist: the distance matrix between all grid locations and the sink
%   params: necessary basic parameters
%
% Return:
%   pwr: the vector of power at all grid locations

function pwr = getPwr(sol, N, dist, params)
% initialization
N_cnt = size(N, 1);     % the number of grid locations

% calculate sensing power
Psen = params.eta * params.Es * sol.s;

% calculate communication power
Pcomm = zeros(N_cnt, 1);
for i = 1:N_cnt
    for j = 1:N_cnt
        if i == j % avoid the same node
            continue;
        end
        fij_idx = (i-1) * N_cnt + j;
        if sol.fij > 0
            % add transmission power to node i
            Pcomm(i) = Pcomm(i) + getPtx(dist(i, j)) * sol.fij(fij_idx) / ...
                params.B;
            % add reception power to node j
            Pcomm(j) = params.Prx * sol.fij(fij_idx) / params.B;
        end
    end
    if sol.fiB(i) > 0
        % add transmission power for node i to sink
        Pcomm(i) = Pcomm(i) + getPtx(dist(i, N_cnt+1)) * sol.fiB(i) / ...
            params.B;
    end
end

% add together
pwr = params.P0 + Psen + Pcomm;
pwr = pwr .* sol.x;        % filter out the power at non-deployed sites
end