%% Reliability-Driven SRIGH algorithm
% Args:
% N: struct of the grid locations
% O: list of targets to monitor
% dist: distance matrix between grid locations and the sink
% params: necessary basic parameters
% rdsrighparams: specific parameters for SRIGH
% 
% Return:
%   sol.fval: optimal value of the objective function
%   sol.exitflag: exit flag showing whether the problem is feasible
%   sol.x: binary vector of optimal node placement
%   sol.s: binary vector optimal sensor placement
%   sol.fij: float vector of flows between grid locations
%   sol.fiB: float vector of flows between grid locations and the sink
function sol = RDSRIGH(N, O, dist, params, rdsrighparams)
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

% Matrix of increased power at node i if link (i, j) is selected
P_inc = zeros(N_cnt+1, N_cnt+1);
for i=1:N_cnt
    for j=i+1:N_cnt+1
        P_inc(i, j) = (getPtx(dist(i, j)) + params.Prx) * params.eta * ...
                params.G / params.B;
        P_inc(j, i) = P_inc(i, j); % P_inc is symmetric
    end
end

% create the graph
G = create_graph(P_inc, N, dist, params, rdsrighparams);
    
% begin selection process
iter = 0;
while sum(T) > 0  
    %% select a sensor and determine the shortest routing path 
    % from that sensor to the sink
    % calculate the current shortest path tree of all nodes to the sink
    [SP, D] = shortestpathtree(G, 1:N_cnt, N_cnt+1);
    
    % calculate the new coverage benefit offered by each node i
    P_bd = horzcat(N(:).Pi);
    new_cover = sum(O_cover' & T, 1) .* P_bd;
    
    % calculate the benefit of placing node at unplaced node i
    benefit = new_cover .* D;
    flag = logical(s == 0);
    benefit = benefit' .* flag;
    
    % get the best benefit value and node index
    [benefit_best, node_best] = max(benefit);
    
    if benefit_best <= 0
        fprintf('There exist targets cannot be covered! Error!\n');
        sol.exitflag = -1;
        return;
    end
    
    % path_best: the list of the nodes in the shortest path
    path_best = shortestpath(G, node_best, N_cnt+1);
    
   %% update sensor placement and unsatisfied coverage 
    % given the new placed sensor
    s(node_best) = 1;
    for j=1:N_o
        if O_cover(node_best, j) > 0
            q(j) = q(j) - 1;
            if q(j) <= 0
                q(j) = 0;
                T(j) = 0;
            end
        end
    end
    
    %% update the graph regarding the nodes on the shortest
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
                    edge_ed = edges(edge_idx, 2);
                    G = rmedge(G, st_idx, edges(edge_idx, 2));
                    new_weight = rdsrighparams.w2 * (P_inc(st_idx, edge_ed) / ...
                        max(1e-10, (N(st_idx).Pi - params.P0 - ...
                        params.Es * params.eta * s(st_idx))));
                    G = addedge(G, st_idx, edges(edge_idx, 2), new_weight);
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
    
    %sol.fval = sum(x);
    %sol.exitflag = 1;
    %sol.x = x;
    %sol.s = s;
    %sol.fij = fij;
    %sol.fiB = fiB;
    %plot_solution(N, O, c, sol, params.S_r, [1000, 1000], sprintf('%d', iter));
    iter = iter + 1;
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
%   P_inc: increased power at i if link (i, j) is selected
%   N: struct of the grid locations
%   dist: distance matrix between grid locations and the sink
%   params: necessary basic parameters
%   rdsrighparams: specific parameters for SRIGH
%
% Return:
%   G: the graph for finding shortest path

function G = create_graph(P_inc, N, dist, params, rdsrighparams)
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
            weight_ij = rdsrighparams.w1 + rdsrighparams.w2 * ...
                (P_inc(i, j)/max(1e-10, (N(i).Pi - params.P0)));
            weight_ji = rdsrighparams.w1 + rdsrighparams.w2 * ...
                (P_inc(j, i)/max(1e-10, (N(j).Pi - params.P0)));
            weights = [weights, weight_ij, weight_ji];
        end
    end
    if dist(i, N_cnt+1) <= params.C_r
        % if connectable, add pair [i, sink] to [st, ed]
        st = [st, i]; ed = [ed, N_cnt+1];
        % update weight from i to sink
        weight_iB = rdsrighparams.w1 + rdsrighparams.w2 * ...
                (P_inc(i, N_cnt+1)/max(1e-10, (N(i).Pi - params.P0)));
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
