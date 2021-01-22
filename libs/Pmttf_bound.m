%% Search for the power bound to satisfy MTTF constraints
%  considering temperature distribution over time.
% Args:
%   rel: the struct of reliability options and targets
%   N: the struct of grid locations
%
% Return:
%   P_mttfi: list of expected upper bound for average power to satisfy MTTF
%            under the current temperature distribution

function P_mttfi = Pmttf_bound(rel, N)
% initialization
N_cnt = size(N, 1);         % number of grid locations
P_mttfi = zeros(N_cnt, 1);  % power bounds imposed by SoH
N_bin = length(N(1).Tcen);  % number of temperature bins
eps = 1e-4;                 % acceptable precision of the output power bound
% start the binary search
for i = 1:N_cnt
    pwr_lb = 0.0;
    pwr_ub = N(i).Ri;
    while pwr_ub - pwr_lb > eps
        pwr_cur = 0.5 * (pwr_lb + pwr_ub);
        MTTF_cur = 0.0;
        % calculate the expectation of SoH over temperature distribution
        for j = 1:N_bin
            Tcellj = amb2core(N(i).Tcen(j), pwr_cur);
            MTTF_cur = MTTF_cur + N(i).Tcnt(j) * mttf(Tcellj);
        end
        %fprintf('MTTF_cur: %f pwr_cur: %f\n', MTTF_cur, pwr_cur);
        % decide the next power lower bound or upper bound
        if MTTF_cur == rel.MTTFref % not really possible
            break;
        elseif MTTF_cur > rel.MTTFref % MTTF is save, current power is low
            pwr_lb = pwr_cur;
        else % MTTF violates the bound, current power is too high
            pwr_ub = pwr_cur;
        end
    end
    %fprintf('Final pwr for %d: %f\n', i, 0.5 * (pwr_lb + pwr_ub));
    %P_mttfi(i) = 0.5 * (pwr_lb + pwr_ub);
    P_mttfi(i) = pwr_ub;
end
end

