% main.m
clc;
clear;
close all;
warning('off','all');
addpath('./lldistkm/');

%% initialize grid map
% pre-process solar and temperature data
fprintf('start pre-processing...\n');
folder = './solardata/';
f_list = dir(append(folder, '*.csv'));
dataT_sav = append(folder, 'dataT.csv');
if exist(dataT_sav, 'file')
    dataT = readtable(dataT_sav);
else
    % if no file exists yet, do pre-process and save it to file
    dataT = preprocess(f_list, dataT_sav);
end

% get data source distribution
bubbleplot_wsize(dataT.lat, dataT.lon, dataT.temp_avg, 'data source');
origin = [min(dataT.lat), min(dataT.lon)];
yScalekm = lldistkm([min(dataT.lat), min(dataT.lon)], [max(dataT.lat), min(dataT.lon)]);
xScalekm = lldistkm([min(dataT.lat), min(dataT.lon)], [min(dataT.lat), max(dataT.lon)]);

% set the size and granularity of the grid space 
xScalem_target = 100;      % m
yScalem_target = 100;      % m
N_x = 3;
N_y = 3;
N_cnt = N_x * N_y;          % number of grid points
Unit_x = floor(xScalem_target / (N_x - 1));
Unit_y = floor(yScalem_target / (N_y - 1));

% the factor to transform from grid (m) to lat and lon
% x - longitude, y - latitude
x_transform = (max(dataT.lon) - min(dataT.lon)) / xScalem_target;
y_transform = (max(dataT.lat) - min(dataT.lat)) / yScalem_target;

% generate the grid candidate set N 
% with their x, y coordinates and temperature and DNI
empty_point.position = [];
empty_point.Ri = [];
empty_point.Ti = [];
N = repmat(empty_point, N_cnt, 1);

for j = 0:N_y-1
    for i = 0:N_x-1
        N(i + j * N_x + 1).position = [Unit_x * i, Unit_y * j];
        x_lon = Unit_x * i * x_transform + origin(2);
        y_lat = Unit_y * j * y_transform + origin(1);
        % find the closest location in dataT
        [minValue, dataT_idx] = min(abs(dataT.lat - y_lat) + ...
            abs(dataT.lon - x_lon));
        % assign the corresponding recharging power in W
        % 0.1 conversion efficiency, 0.01m2 solar panel, w/m2 radiation
        N(i + j * N_x + 1).Ri = 0.1 * 0.01 * dataT.dni_avg(dataT_idx);
        % assign the corresponding temperature in Celsius
        N(i + j * N_x + 1).Ti = dataT.temp_avg(dataT_idx);
    end
end

%% prepare for constraints
S_r = 40;           % sensing range in m
C_r = 50;           % communication range in m
N_o = 2;            % number of PoIs
K = 1;              % K-coverage
G = 100;            % generated bytes for each sample
B = 2000;           % bandwidth B/s
L_eta = 0.001;       % lower bound for duty cycle
% randomly generate PoIs to monitor
O = repmat([], N_o, 2);
for i = 1:N_o
    O(i, 1) = unifrnd(0, xScalem_target);
    O(i, 2) = unifrnd(0, yScalem_target);
end
% location of the sink
c = [xScalem_target/2, yScalem_target/2];

% the binary vectors to locate each type of variable
v_cnt = N_cnt + N_cnt + N_cnt + N_cnt^2 + N_cnt; % x, s, eta, fij, fiB
% useful parameters
% cnt: number of variables in this type.
% base: number of variables before this type.
x_cnt = N_cnt; x_base = 0;
s_cnt = N_cnt; s_base = N_cnt;
eta_cnt = N_cnt; eta_base = 2*N_cnt;
fij_cnt = N_cnt^2; fij_base = 3*N_cnt;
fiB_cnt = N_cnt; fiB_base = 3*N_cnt + N_cnt^2;

%% call the solver
% objective
v_x = [ones(1, x_cnt), zeros(1, v_cnt-s_base)];
fun = @(x) v_x*x;    % Objective Function f(x)

% linear inequality coverage constraints padded to N_o * v_cnt
A = -cover_matrix(O, N, S_r, v_cnt);
b = -repmat(K, N_o, 1);
% linear equality constraint
% set s to zeros
Aeq_set = [zeros(s_cnt, x_cnt), eye(x_cnt), zeros(x_cnt, v_cnt-eta_base)];
beq_set = zeros(s_cnt, 1);
% infeasible flows
[Aeq_infeasible, beq_infeasible] = flow_infeasible(N, v_cnt, fij_base, C_r);
Aeq = [Aeq_set; Aeq_infeasible];
beq = [beq_set; beq_infeasible];

% nonlinear constraints
nlcon = @(x) [getPower(x, G, B, N, c); flow(x, G, N_cnt)];
nlrhs = [vertcat(N(:).Ri); zeros(N_cnt+1, 1)];
% -1 for <=, 0 for ==, +1 >=  
nle = [repmat(-1, N_cnt, 1); zeros(N_cnt+1, 1)];

% bounds
lb = [zeros(eta_base, 1); repmat(L_eta, eta_cnt, 1); zeros(v_cnt-fij_base, 1)];
ub = [ones(fij_base, 1); Inf(v_cnt-fij_base, 1)];

% integer constraints
nB = x_cnt + s_cnt;                 % Number of Binary Variables
nC = eta_cnt + fij_cnt + fiB_cnt;   % Number of Continuous Variables
nI = 0;                             % Number of Integer Variables

% build xtype vector
xtype = [repmat('B', 1, nB), repmat('C', 1, nC), repmat('I', 1, nI)];

% initial guess
x0 = zeros(v_cnt, 1);

% create OPTI Object
Opt = opti('fun', fun, 'nlmix', nlcon, nlrhs, nle, 'ineq', A, b, ...
        'bounds', lb, ub, 'eq', Aeq, beq, 'xtype', xtype)

% solve the MINLP problem
[x,fval,exitflag,info] = solve(Opt,x0)
% extract flow matrix fij (N_cnt * N_cnt) from x
T_fij = [zeros(N_cnt^2, 3*N_cnt), eye(N_cnt^2), zeros(N_cnt^2, N_cnt)];
fij = T_fij * x;
fij = reshape(fij, N_cnt, N_cnt)'; % (i, j) is flow from i to j 
fij

% plot the solution
place = x(1:N_cnt); % first N_cnt binary variables is the placement solution
plot_solution(N, O, place, S_r);

%% plot functions
% plot the locations
function bubbleplot_wsize(lat, lon, sizedata, title)
    figure;
    geobubble(lat, lon, sizedata, 'Title', title);
    ax = gca; % get current axes
    ax.FontSize = 16;
    %geobasemap streets-light; % set base map style
end

% plot the solution in the grid space
function plot_solution(N, O, x, S_r)
    color_map = [[0 0.4470 0.7410]; ...       % blue
        [0.8500 0.3250 0.0980]; ...           % orange
        [0.4660 0.6740 0.1880]];              % green
    nodes = vertcat(N.position);              % grid points
    nodes = vertcat(nodes, O);                % PoIs
    sz_nodes = repmat(10, size(nodes, 1), 1); % const size for nodes
    color_idx = vertcat(x, repmat(2, size(O, 1), 1)) + 1;
    color_nodes = vertcat(color_map(color_idx, :));
    figure;
    scatter(nodes(:, 1), nodes(:, 2), sz_nodes, color_nodes);
end