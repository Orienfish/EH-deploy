%% Temperature-Aware Two-Stage Heuristic
% Args:
%   N: struct of the grid locations
%   O: list of targets to monitor
%   dist: distance matrix between grid locations and the sink
%   params: necessary basic parameters
%   tatshparams: specific parameters for TATSH
%
% Return:
%   sol.fval: optimal value of the objective function
%   sol.exitflag: exit flag showing whether the problem is feasible
%   sol.x: binary vector of optimal node placement
%   sol.s: binary vector optimal sensor placement
%   sol.fij: float vector of flows between grid locations
%   sol.fiB: float vector of flows between grid locations and the sink
%   sol.Pi: left-hand side and right-hand side of power inequalities
%   sol.cov: coverage of each target

function sol = TATSH(N, O, dist, params, tatshparams)
% initialization
N_cnt = size(N, 1);             % number of grid locations
N_o = size(O, 1);               % number of targets to monitor
%x = zeros(N_cnt, 1);            % binary vector of node placement
s = zeros(N_cnt, 1);            % binary vector of sensor placement
q = repmat(params.K, N_o, 1);   % coverage requirement of each target
T = ones(N_o, 1);               % binary vector of uncovered targets

%% stage 1: select sensor nodes
while sum(T) > 0 % loop continues if there is uncovered target
    % init the weight vector in selecting sensor nodes
    w = -Inf(N_cnt, 1);
    % calculate the weight of each unselected sensor
    for i=1:N_cnt
        if s(i) == 0
            O_cover = cover_targets(N(i).position, O, params.S_r);
            %disp(O_cover');
            %disp(T');
            w(i) = sum(O_cover & T) / N(i).Ti;
            %fprintf('%d: %f\n', i, w(i));
        end
    end
    %disp(w');
    % select the sensor with maximum weight in this round
    [maxval, maxidx] = max(w);
    if maxval <= 0 % exist targets cannot be covered, return failure
        sol.exitflag = -1;
        return;
    end
    % update q and T with selection maxidx
    O_cover = cover_targets(N(maxidx).position, O, params.S_r);
    new_cover = O_cover & T;
    for j=1:N_o
        if new_cover(j) > 0
            q(j) = q(j) - 1;
            if q(j) <= 0
                T(j) = 0; 
            end
        end
    end
    s(maxidx) = 1;   
end
x = s;


%% stage 2: select relay nodes
% update the solution value
sol.fval = sum(x);
sol.exitflag = 1;
sol.x = x;
sol.s = s;
sol.fij = zeros(N_cnt^2, 1);
sol.fiB = zeros(N_cnt, 1);
end

%% get the covered target by node i
function O_cover = cover_targets(pos_i, O, S_r)
    N_o = size(O, 1);               % number of targets to monitor
    O_cover = zeros(N_o, 1);
    for j=1:N_o
        if norm(O(j, :) - pos_i) <= S_r
            O_cover(j) = 1;
        end
    end
end