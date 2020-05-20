%% Two-Stage Heuristic from the directional paper
% Args:
%   N: struct of the grid locations
%   O: list of targets to monitor
%   dist: distance matrix between grid locations and the sink
%   params: necessary basic parameters
%   tshparams: specific parameters for TSH
%
% Return:
%   sol.fval: optimal value of the objective function
%   sol.exitflag: exit flag showing whether the problem is feasible
%   sol.x: binary vector of optimal node placement
%   sol.s: binary vector optimal sensor placement
%   sol.fij: float vector of flows between grid locations
%   sol.fiB: float vector of flows between grid locations and the sink

function sol = TSH(N, O, dist, params, tshparams)
% initialization
N_cnt = size(N, 1);             % number of grid locations
N_o = size(O, 1);               % number of targets to monitor
%x = zeros(N_cnt, 1);            % binary vector of node placement
s = zeros(N_cnt, 1);            % binary vector of sensor placement
q = repmat(params.K, N_o, 1);   % coverage requirement of each target
T = ones(N_o, 1);               % binary vector of uncovered targets
fij = zeros(N_cnt^2, 1);        % flow matrix between grids
fiB = zeros(N_cnt, 1);          % flow vector to the sink

% the cost vector of adding a solar panel
cost = tshparams.w2 * (getPtx(params.C_r) + params.Prx) * params.eta * ...
    params.G / params.B ./ vertcat(N(:).Ri);

%% stage 1: select sensor nodes
while sum(T) > 0 % loop continues if there is uncovered target
    % init the weight vector in selecting sensor nodes
    w = -Inf(N_cnt, 1);
    % calculate the weight of each unselected sensor
    for i=1:N_cnt
        if s(i) == 0
            O_cover = cover_targets(N(i).position, O, params.S_r);
            w(i) = sum(O_cover & T);
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

%% stage 2: select relay nodes for the placed sensor one by one
% create the directed network graph
G = create_graph(x, cost, N, dist, params, tshparams);

% find the shortest path from all selected sensors to the sink
senidx = find(s > 0); % get the indexes of selected sensors
SP = shortestpathtree(G, senidx, N_cnt+1);
% extract the [st, ed] nodes pair in the shortest path
pair = reshape(SP.Edges.EndNodes(:), [], 2);
for i=1:size(pair, 1)
    st = pair(i, 1); ed = pair(i, 2);
    x(st) = 1; % update relay node placement
    % update flow matrix
    if ed <= N_cnt % sending to another grid
        fij_idx = (st-1) * N_cnt + ed;
        fij(fij_idx) = fij(fij_idx) + params.eta * params.G;
    else % sending to the sink
        fiB(st) = fiB(st) + params.eta * params.G;
    end
end

%% update the solution value
sol.fval = sum(x);
sol.exitflag = 1;
sol.x = x;
sol.s = s;
sol.fij = fij;
sol.fiB = fiB;
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

%% create the network graph
% Args:
%   x: binary vector of current node placement
%   cost: cost vector of adding a solar panel
%   N: struct of the grid locations
%   dist: distance matrix between grid locations and the sink
%   params: necessary basic parameters
%   tshparams: specific parameters for TSH
%
% Return:
%   G: the graph for finding shortest path

function G = create_graph(x, cost, N, dist, params, tshparams)
N_cnt = size(N, 1);
st = [];
ed = [];
weights = [];
for i=1:N_cnt
    for j=i+1:N_cnt
        if dist(i, j) <= params.C_r
            % if connectable, add pair [i, j] and [j, i] to [st, ed]
            st = [st, i, j]; ed = [ed, j, i];
            % update weight from i to j and from j to i
            weight_ij = tshparams.w1 * (x(i) == 0) + tshparams.w2 * cost(i);
            weight_ji = tshparams.w1 * (x(j) == 0) + tshparams.w2 * cost(j);
            weights = [weights, weight_ij, weight_ji];
        end
    end
    if dist(i, N_cnt+1) <= params.C_r
        % if connectable, add pair [i, sink] to [st, ed]
        st = [st, i]; ed = [ed, N_cnt+1];
        % update weight from i to sink
        weight_iB = tshparams.w1 * (x(i) == 1) + tshparams.w2 * cost(i);
        weights = [weights, weight_iB];
    end
end
G = digraph(st, ed, weights);
%G.Edges
end
