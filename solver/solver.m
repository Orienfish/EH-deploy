%% Call the solver the calculate the optimal solution 
% Args:
%   N: struct of the grid locations
%   O: list of targets to monitor
%   dist: distance matrix between grid locations and the sink
%   params: necessary basic parameters
%   rel: the reliability options and targets
%
% Return:
%   sol.fval: optimal value of the objective function
%   sol.x: binary vector of optimal node placement
%   sol.s: binary vector optimal sensor placement
%   sol.fij: float vector of flows between grid locations
%   sol.fiB: float vector of flows between grid locations and the sink
%   sol.Pi: float vector of power in W for each node
%   sol.cov: coverage of each target

function sol = solver(N, O, dist, params, rel)

addpath('./libs');
%% Prepare the constraints for the solver
fprintf('Preparing the constraints...\n');

% the binary vectors to locate each type of variable
N_cnt = size(N, 1); % number of grid locations 
v_cnt = N_cnt + N_cnt + N_cnt^2 + N_cnt; % x, s, fij, fiB
% parameters of the format of variables
% cnt: number of variables in this type.
% base: number of variables before this type.
xform.x_cnt = N_cnt; xform.x_base = 0; 
xform.x_end = xform.x_base + xform.x_cnt;
xform.s_cnt = N_cnt; xform.s_base = N_cnt; 
xform.s_end = xform.s_base + xform.s_cnt;
xform.fij_cnt = N_cnt^2; xform.fij_base = 2*N_cnt;
xform.fij_end = xform.fij_base + xform.fij_cnt;
xform.fiB_cnt = N_cnt; xform.fiB_base = 2*N_cnt + N_cnt^2;
xform.fiB_end = xform.fiB_base + xform.fiB_cnt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective Function f*x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = [ones(1, xform.x_cnt), zeros(1, v_cnt-xform.x_end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear equality constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flow conservation constraints
[Afc, bfc] = flow_conserve(params.G, params.eta, v_cnt, xform, dist, params.C_r);
% comple connectivity constraints
v_fiB = [zeros(1, xform.fiB_base), ones(1, xform.fiB_cnt)];
v_s = [zeros(1, xform.x_cnt), ones(1, xform.s_cnt), ...
    zeros(1, v_cnt - xform.s_end)];
Acc = v_fiB - params.eta * params.G * v_s; % 1*1
bcc = 0;
% stack all equality constraints
Aeq = [Afc; Acc];
beq = [bfc; bcc];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear inequality constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coverage constraints
[Ac, bc] = coverage(O, N, params.S_r, v_cnt, params.K);
% xi >= si
% get the transformation matrix to extract x and s
T_x = [eye(xform.x_cnt), zeros(xform.x_cnt, v_cnt-xform.x_end)];
T_s = [zeros(xform.s_cnt, xform.x_cnt), eye(xform.s_cnt), ...
    zeros(xform.s_cnt, v_cnt-xform.s_end)];
Axs = T_s - T_x;
bxs = zeros(N_cnt, 1);
% sigma_j fij <= maxf * xi
[Abound, bbound] = flow_bound(v_cnt, xform, dist, params.C_r, params.maxf);
% power inequality
[Apwr, bpwr] = power_eno(params.eta, params.B, v_cnt, xform, N, dist, ...
    params.C_r, rel);
% stack all inequality constraints
Aineq = [Ac; Axs; Abound; Apwr];
bineq = [bc; bxs; bbound; bpwr];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lb = zeros(v_cnt, 1);
ub = [ones(xform.fij_base, 1); repmat(params.maxf, v_cnt-xform.fij_base, 1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% integer constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nB = xform.x_cnt + xform.s_cnt;           % Number of Binary Variables
nC = xform.fij_cnt + xform.fiB_cnt;       % Number of Continuous Variables
nI = 0;                                   % Number of Integer Variables

% build ctype vector
ctype = [repmat('B', 1, nB), repmat('C', 1, nC), repmat('I', 1, nI)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial guess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = zeros(v_cnt, 1);

%% call the solver
fprintf('Calling CPLEX...\n');
options = cplexoptimset;
options.Display = 'Off';
[x, fval, exitflag, output] = cplexmilp(f, Aineq, bineq, Aeq, beq, ...
    [], [], [], lb, ub, ctype, x0);
output

% fill the solution
sol.fval = fval;
sol.x = x(xform.x_base+1: xform.x_end);
sol.s = x(xform.s_base+1: xform.s_end);
sol.fij = x(xform.fij_base+1: xform.fij_end);
sol.fiB = x(xform.fiB_base+1: xform.fiB_end);
global P0;
sol.Pi = Apwr * x + repmat(P0, N_cnt, 1);
sol.cov = Ac * x;

end