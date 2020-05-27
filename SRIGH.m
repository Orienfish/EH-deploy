%% SRIGH algorithm
% Cparams: basic parameters
% sparams: parameters for srigh
function [res] = SRIGH(Cparams, sparams)
% Cparams.M: number of deployment sites
% Cparams.N: number of targets
% Cparams.S_r: sensing radius
% Cparams.C_r: communication range in m
% Cparams.K: coverage requirements
% Cparams.dist: distance matrix between grid locations and the sink
% Cparams.Prx: reception power (W)
% Cparams.eta: sampling frequency
% Cparams.G: generated bytes for each sample
% Cparams.B: bandwidth B/s
% 
% sparams.A: surface area of the solar panel
% sparams.x: average solar irradiance
% sparams.O: deployment sites, struct
% sparams.T: targets
%
% res.exitflag: exit flag showing whether the problem is feasible
% res.x: binary vector of optimal node placement
% res.s: binary vector optimal sensor placement
% res.fij: float vector of flows between grid locations
% res.fiB: float vector of flows between grid locations and the si

s = zeros(Cparams.M, 1);            % binary vector of sensor placement
x = zeros(Cparams.M, 1);            % binary vector of node placement
q = repmat(Cparams.K, Cparams.N, 1);  % coverage requirement of each target
fij = zeros((Cparams.M)^2, 1);        % flow matrix between grids
fiB = zeros(Cparams.M, 1);          % flow vector to the sink
% eff_s = get_efficient_sensors(sparams.O, sparams.T, Cparams.S_r, ...,
%                               Cparams.M, Cparams.N);
% gamma = eff_s;   % binary vector of sensors that could cover some targets
%                  % and no sensors have been placed

% T' <- the set for targets
% 0 if the covering requirement has been fulfilled, 1 if not
T_prime = ones(1, Cparams.N);

% qi' <- qi
% targets with their unfulfilled covering requirements
qi_prime = q;

% omega_ik <- binary matrix storing whether the sensor in row i could
% cover the target in column j
o_cover = [];
O = sparams.O;
% gamma: binary vector of sensors that could cover some targets
gamma = zeros(Cparams.M, 1);
for i = 1:Cparams.M
    o_ik = cover_targets(O(i).position, sparams.T, Cparams.S_r);
    o_cover = vertcat(o_cover, o_ik');
    gamma(i) = logical(sum(o_ik));
end

% the cost vector of adding a solar panel
cost = sparams.w2 * (getPtx(Cparams.C_r) + Cparams.Prx) * Cparams.eta * ...
    Cparams.G / Cparams.B ./ vertcat(O(:).Ri);

% create the graph
G = create_graph(cost, sparams.O, Cparams.dist, Cparams.C_r, ...
        sparams.w1);
    
% plot(G,'EdgeLabel',G.Edges.Weight)
% begin selection process
while sum(T_prime) > 0
    % check whether to exit loop
    for i = 1:Cparams.N
        % If there's a column have sum 0, it means that none of the sensors
        % could satisfy covering this target.
        % If this target is one of the targets whose covering requirement
        % hasn't been met, then it's not possible to cover this target.
        if(sum(o_cover(:,i)) == 0 && T_prime(i) == 1)
            res.exitflag = -1;
            return
        end
    end
    
    %% stage 1: selecting a target
    % h <- argmin(j in T') {(sum ((i in gamma) {omega i, j}))/qj }
    minj = 1; % index of target with min rho_i/q'_i
    minval = Inf; % rho_i/q'_i of that target
    for j = 1:Cparams.N
        % skip if the target has fulfilled the coverage requirement
        if T_prime(j) == 0
            continue
        end
        num_sensors = 0;
        % sum (i in gamma) {omega(i, target)}
        for i = 1:Cparams.M
            if gamma(i) == 1
                num_sensors = num_sensors + o_cover(i, j);
            end
        end
        % rho_i/q'_i of current target
        temp = num_sensors / qi_prime(j);
        % update target
        if temp < minval
            minval = temp;
            minj = j;
        end
    end
    h = minj;
    
    %% stage 2: select a sensor to cover the target
    g_s = 0;  % best benefit
    i_s = -1; % sensor with best benefit
    P_s = []; % path with best benefit
    
    % get best benefit, sensor with best benefit, and the path to achieve
    % the benefit
    for i = 1:Cparams.M
        % skip the sensor if it cannot cover the selected target
        if o_cover(i, h) == 0 || gamma(i) == 0
            continue
        end
        
        % P: the list of the nodes in the shortest path
        % d: the length of the shortest path
        [P, d] = shortestpath(G, i, Cparams.M+1);
        
        % g <- |Ti intersect T'|/(d + cost(i))
        % first calculate numerator
        % Ti: binary vector of targets that could be covered by current sensor
        Ti = o_cover(i, :);
        intersection = Ti & T_prime;
        numerator = sum(intersection);

        % then calculate denominator
        denominator = d + cost(i);
        % If i (current sensor) is already selected, then don't need to add
        % omega1 when calculating the cost.
        denominator = denominator + sparams.w1 * (x(i) == 0);
        
        % calculate benefit
        g = numerator / denominator;
        % update best benefit, sensor with best benefit, and thae
        % corresponding path
        if g > g_s
            g_s = g;
            i_s = i;
            P_s = P;
        end
    end
    
    %% stage 3: update target set and undeployed sensor set(gamma)
    % build the set of targets who could be covered by the selected sensor
    % T_is: binary vector, 1 if the target could be covered, 0 if not
    T_is = o_cover(i_s, :);
    % intersection: unfulfilled targets that could be covered by the selected sensor
    intersection = T_is & T_prime;
    
    % update the unfulfilled target set
    for i = 1:Cparams.N
        % skip if the target cannot be covered
        if intersection(i) == 0
            continue
        end
        % decrement the coverage requirement
        qi_prime(i) = qi_prime(i) - 1;
        
        % update T_prime (set of unfulfilled targets) if the coverage
        % requirement has been fulfilled
        if qi_prime(i) <= 0
            T_prime(i) = 0;
        end
    end
    
    % update undeployed sensor set
    gamma(i_s) = 0;
    
    % update selected sensor sets
    s(i_s) = 1;
    %x(i_s) = 1;
    %% stage 4: update the graph
    edges = G.Edges.EndNodes; % M by 2 matrix
    endpoints = edges(:, 2); % get the second column, the endpoints
    
    % P_s contains the sequence of nodes
    % traversing from 1st to second last nodes in P_s
    % (O_i -> O_j) == arc <O_i, O_j>
    for i = 1:(size(P_s, 2)-1)
        O_i = P_s(i);
        O_j = P_s(i + 1);
        
        % add the newly-selected node in Fij/FiB
        if O_j <= Cparams.M % sending to another grid
            fij_idx = (O_i - 1) * Cparams.M + O_j;
            fij(fij_idx) = fij(fij_idx) + Cparams.eta * Cparams.G;
        else % sending to the sink
            fiB(O_i) = fiB(O_i) + Cparams.eta * Cparams.G;
        end
        
        % if O_i is not in the "selected node" set, then need to add it in
        % and update the graph weight
        if x(O_i) == 0
            x(O_i) = 1;
            
            % update weight
            for j = 1:size(endpoints, 1)
                % endpoint is O_i
                if (endpoints(j) == O_i)
                    % get the start point
                    % arc: <startpoint, O_i>
                    startpoint = edges(j, 1);
                    weight = cost(O_i);
                    % update the edge
                    G = rmedge(G, startpoint, O_i);
                    G = addedge(G, startpoint, O_i, weight);
                end
            end
        end
    end
end

res.fval = sum(x);
res.s = s;
res.fij = fij;
res.x = x;
res.fiB = fiB;
res.exitflag = 1;
end
%% create the network graph
% Args:
%   cost: cost vector of adding a solar panel
%   N: struct of the grid locations
%   dist: distance matrix between grid locations and the sink
%   C_r: communication radius
%   w1: cost for adding a new node
%   w2: cost for adding per area of solar panel
%
% Return:
%   G: the graph for finding shortest path

function G = create_graph(cost, N, dist, C_r, w1)
N_cnt = size(N, 1);
st = [];
ed = [];
weights = [];
for i=1:N_cnt
    for j=i+1:N_cnt
        if dist(i, j) <= C_r
            % if connectable, add pair [i, j] and [j, i] to [st, ed]
            st = [st, i, j]; ed = [ed, j, i];
            % update weight from i to j and from j to i
            weight_ij = w1 + cost(j);
            weight_ji = w1 + cost(i);
            weights = [weights, weight_ij, weight_ji];
        end
    end
    if dist(i, N_cnt+1) <= C_r
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
