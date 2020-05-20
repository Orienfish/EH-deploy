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
            w(i) = sum(O_cover & T) / N(i).Ti;
        end
    end
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

% calculate the current power Pi for all grid locations
P0 = 0.01;                      % sleep power (W)
Es = 0.2*0.2;                   % sensing power (W*s)
Prx = 0.1;                      % reception power (W)
P_cur = P0 + s * Es * params.eta;  % list of current power - only sensing
%disp(Pi);

%% stage 2: select relay nodes for the placed sensor one by one
% create the directed network graph
st = [];
ed = [];
weights = [];
for i=1:N_cnt
    for j=i+1:N_cnt
        if dist(i, j) <= params.C_r
            % if connectable, add pair [i, j] and [j, i] to [st, ed]
            st = [st, i, j]; ed = [ed, j, i];
            % increased transmission power due to placing relay node at i or j
            P_inc = (getPtx(dist(i, j)) + Prx) * params.eta * params.G / params.B;
            weight_ij = tatshparams.w1 * (s(i) == 0) + ...
                tatshparams.w2 * (1/max(1e-10, (N(i).Ri - P_cur(i) - P_inc)));
            weight_ji = tatshparams.w1 * (s(j) == 0) + ...
                tatshparams.w2 * (1/max(1e-10, (N(j).Ri - P_cur(j) - P_inc)));
            weights = [weights, weight_ij, weight_ji];
        end
    end
    if dist(i, N_cnt+1) <= params.C_r
        % if connectable, add pair [i, sink] to [st, ed]
        st = [st, i]; ed = [ed, N_cnt+1];
        % increased transmission power due to placing relay node at i
        P_inc = (getPtx(dist(i, N_cnt+1)) + Prx) * params.eta * params.G / params.B;
        weight_iB = tatshparams.w1 * (s(i) == 1) + ...
            tatshparams.w2 * (1/max(1e10, (N(i).Ri - P_cur(i) - P_inc)));
        weights = [weights, weight_iB];
    end
end
G = digraph(st, ed, weights);
disp(G);

% start to choose the shortest path to sink for each sensor in stage 1
prior = dist(:, N_cnt+1) .* s; % the priority of each placed sensor
prior(prior == 0) = Inf;       % set the zero item in prior to Inf
                               % so we can start from the node with min
                               % prior value
%disp(prior);


%% update the solution value
sol.fval = sum(x);
sol.exitflag = 1;
sol.x = x;
sol.s = s;
sol.fij = zeros(N_cnt^2, 1);
sol.fiB = zeros(N_cnt, 1);
end

%% get the covered target by node i
% Args:
%   pos_i: the location of the ith node
%   O: list of targets to monitor
%   S_r: sensing radius
%
% Return:
%   O_cover: binary vector showing which targets can be covered by i

function O_cover = cover_targets(pos_i, O, S_r)
    N_o = size(O, 1);               % number of targets to monitor
    O_cover = zeros(N_o, 1);
    for j=1:N_o
        if norm(O(j, :) - pos_i) <= S_r
            O_cover(j) = 1;
        end
    end
end

%% Get the transmission power from distance
% Args:
%   dist_m: the distance in m
%
% Return:
%   Ptx: the transmission power in W

function Ptx = getPtx(dist_m)
Pto = 0.22;                 % transmission power baseline
alpha = 3.5;                % distance coefficient
beta = 1e-7;                % distance coefficient

Ptx = Pto + beta * dist_m ^ alpha;
end
