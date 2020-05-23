%% Search for the power bound to satisfy SoH constraints
%  considering temperature distribution over time.
% Args:
%   rel: the struct of reliability options and targets
%   N: the struct of grid locations
%
% Return:
%   P_sohi: list of expected upper bound for average power to satisfy SoH
%           under the current temperature distribution

function P_sohi = Psoh_bound(rel, N)
% initialization
N_cnt = size(N, 1);         % number of grid locations
P_sohi = zeros(N_cnt, 1);   % power bounds imposed by SoH
N_bin = length(N(1).Tcen);  % number of temperature bins
eps = 1e-4;                 % acceptable precision of the output power bound
% start the binary search
for i = 1:N_cnt
    pwr_lb = 0.0;
    pwr_ub = N(i).Ri;
    while pwr_ub - pwr_lb > eps
        pwr_cur = 0.5 * (pwr_lb + pwr_ub);
        SoH_cur = 0.0;
        % calculate the expectation of SoH over temperature distribution
        for j = 1:N_bin
            Tcellj = amb2core(N(i).Tcen(j), pwr_cur);
            SoH_cur = SoH_cur + N(i).Tcnt(j) * soh(Tcellj, rel.T);
        end
        %fprintf('SoH_cur: %f pwr_cur: %f\n', SoH_cur, pwr_cur);
        % decide the next power lower bound or upper bound
        if SoH_cur == rel.SoHref % not really possible
            break;
        elseif SoH_cur > rel.SoHref % SoH is save, current power is low
            pwr_lb = pwr_cur;
        else % SoH violates the bound, current power is too high
            pwr_ub = pwr_cur;
        end
    end
    %fprintf('Final pwr for %d: %f\n', i, 0.5 * (pwr_lb + pwr_ub));
    P_sohi(i) = 0.5 * (pwr_lb + pwr_ub);
end
end

