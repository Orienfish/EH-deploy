%% SRIGH algorithm
% Args:
% N: struct of the grid locations
% O: list of targets to monitor
% dist: distance matrix between grid locations and the sink
% params: necessary basic parameters
% srighparams: specific parameters for SRIGH
% 
% Return:
%   sol.fval: optimal value of the objective function
%   sol.exitflag: exit flag showing whether the problem is feasible
%   sol.x: binary vector of optimal node placement
%   sol.s: binary vector optimal sensor placement
%   sol.fij: float vector of flows between grid locations
%   sol.fiB: float vector of flows between grid locations and the sink
function sol = SRIGH(N, O, dist, params, srighparams)
% initialization
N_cnt = size(N, 1);             % number of grid locations
N_o = size(O, 1);               % number of targets to monitor
s = zeros(N_cnt, 1);            % binary vector of sensor placement
x = zeros(N_cnt, 1);            % binary vector of node placement
q = repmat(params.K, N_o, 1);   % unsatisfied coverage requirement of each target
T = ones(N_o, 1);               % binary vector of uncovered targets
fij = zeros(N_cnt^2, 1);        % flow matrix between grids
fiB = zeros(N_cnt, 1);          % flow vector to the sink

% O_cover: binary matrix indicating whether the i th sensor could
% cover the j th target
O_cover = cover_targets(N, O, params.S_r);

% the cost vector of adding a solar panel
cost = (getPtx(params.C_r) + params.Prx) * params.eta * ...
    params.G / params.B ./ vertcat(N(:).Ri);

% create the graph
G = create_graph(cost, N, dist, params, srighparams);
    
% begin selection process
while sum(T) > 0
    % check whether to exit loop
    % if no more new coverage can be made, then exit
    unplaced_node_idx = logical(s == 0);
    uncover_target_idx = logical(T == 1);
    if sum(O_cover(unplaced_node_idx, uncover_target_idx)) == 0
        fprintf('There exist targets cannot be covered! Error!\n')
        sol.exitflag = -1;
        return;
    end
    
    %% stage 1: selecting a target
    % calculate the number of sensors that can cover each target
    rho = sum(O_cover, 1)';
    % calculate the priority value
    priority = rho ./ q;
    % the target with less number of sensors can cover and more unsatisfied
    % coverage level is given higher priority
    [priority_val, tgt_idx] = min(priority);
    
    %% stage 2: select a sensor to cover the target and determine
    % the shortest routing path from that sensor to the sink
    benefit_best = 0;   % best benefit
    
    % get best benefit, sensor with best benefit, and the path to achieve
    % the benefit
    for i=1:N_cnt
        % skip the sensor if it cannot cover the selected target or
        % there is a sensor already placed at this location
        if O_cover(i, tgt_idx) == 0 || s(i) == 1
            continue
        end
        
        % P: the list of the nodes in the shortest path
        % d: the length of the shortest path
        [P, d] = shortestpath(G, i, N_cnt+1);
        
        % calculate benefit of placing node at i
        % benefit <- |Ti intersect T'|/(d + w1 * (whether a sensor is placed))
        % first calculate numerator: new coverage that has not been covered
        % in previous selections
        new_cover = sum(O_cover(i, :)' & T);
        
        % then calculate the denominator: cost of shortest path plus 
        % new node placement cost
        % If no node has been placed at i, a new-node cost is added
        new_cost = d + srighparams.w1 * (x(i) == 0);
        
        benefit = new_cover / new_cost;
        
        % update best benefit
        if benefit > benefit_best
            benefit_best = benefit;
            node_best = i;
            path_best = P;
        end
    end
    
    %% stage 3: update sensor placement and unsatisfied coverage 
    % given the new placed sensor
    s(node_best) = 1;
    for j=1:N_o
        if O_cover(node_best, j) > 0
            q(j) = q(j) - 1;
            if q(j) <= 0
                q(j) = 1e-4; % avoid q to be zero since it appears as 
                             % denominator in target selection
                T(j) = 0;
            end
        end
    end
    
    %% stage 4: update the graph since the nodes on the shortest
    % path have already been placed
    % the first term in computing weight(i, j) should be removed
    edges = G.Edges.EndNodes; % M by 2 matrix
    
    % update along the shortest path
    for i = 1:(size(path_best, 2)-1)
        st_idx = path_best(i);         % start node on this edge
        ed_idx = path_best(i + 1);     % end node on this edge
        
        % if source node is not in the "selected node" set
        % then add it in and update the graph weight
        if x(st_idx) == 0
            x(st_idx) = 1;       
            % update weight
            for edge_idx = 1:size(edges, 1)
                % if the start point of the current edge is st_idx
                if (edges(edge_idx, 1) == st_idx)
                    % update the edge
                    G = rmedge(G, st_idx, edges(edge_idx, 2));
                    G = addedge(G, st_idx, edges(edge_idx, 2), cost(st_idx));
                end
            end
        end
        
        % update flow matrix
        if ed_idx <= N_cnt % sending to another grid
            fij_idx = (st_idx-1) * N_cnt + ed_idx;
            fij(fij_idx) = fij(fij_idx) + params.eta * params.G;
        else % sending to the sink
            fiB(st_idx) = fiB(st_idx) + params.eta * params.G;
        end
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


%% create the network graph
% Args:
%   cost: cost vector of adding a solar panel
%   N: struct of the grid locations
%   dist: distance matrix between grid locations and the sink
%   params: necessary basic parameters
%   srighparams: specific parameters for SRIGH
%
% Return:
%   G: the graph for finding shortest path

function G = create_graph(cost, N, dist, params, srighparams)
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
            weight_ij = srighparams.w1 + srighparams.w2 * cost(i);
            weight_ji = srighparams.w1 + srighparams.w2 * cost(j);
            weights = [weights, weight_ij, weight_ji];
        end
    end
    if dist(i, N_cnt+1) <= params.C_r
        % if connectable, add pair [i, sink] to [st, ed]
        st = [st, i]; ed = [ed, N_cnt+1];
        % update weight from i to sink
        %weight_iB = w1 + cost(i);
        weight_iB = 0;
        weights = [weights, weight_iB];
    end
end
G = digraph(st, ed, weights);
%G.Edges
end

%% get the covered target by node i
% Args:
%   N: struct of the grid locations
%   O: list of targets to monitor
%   S_r: sensing radius
%
% Return:
%   O_cover: binary matrix with the (i, j)th element showing whether
%            sensor node i can cover targets j

function O_cover = cover_targets(N, O, S_r)
N_cnt = size(N, 1);         % number of potential sensor nodes
N_o = size(O, 1);               % number of targets to monitor
O_cover = zeros(N_cnt, N_o);
for i=1:N_cnt
    for j=1:N_o
        if norm(N(i).position - O(j, :)) <= S_r
            O_cover(i, j) = 1;
        end
    end
end
end
