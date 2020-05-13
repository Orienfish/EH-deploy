% main.m
clc;
clear;
close all;
warning('off','all');
addpath ('/Applications/CPLEX_Studio1210/cplex/matlab/x86-64_osx');
addpath('./lldistkm');

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
N_x = 4;
N_y = 4;
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
G = 32;             % generated bytes for each sample
B = 2000;           % bandwidth B/s
eta = 0.01;         % sampling frequency
% randomly generate PoIs to monitor
O = repmat([], N_o, 2);
for i = 1:N_o
    O(i, 1) = unifrnd(0, xScalem_target);
    O(i, 2) = unifrnd(0, yScalem_target);
end
% location of the sink
c = [xScalem_target/2, yScalem_target/2];

% the binary vectors to locate each type of variable
x_cnt = N_cnt + N_cnt + N_cnt^2 + N_cnt; % x, s, fij, fiB
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
xform.fiB_end = xform.fiB_base + xform.fij_cnt;

%% call the solver
 % get the transformation matrix to extract x and s
T_x = [eye(xform.x_cnt), zeros(xform.x_cnt, x_cnt-xform.x_end)];
T_s = [zeros(xform.s_cnt, xform.x_cnt), eye(xform.s_cnt), ...
    zeros(xform.s_cnt, x_cnt-xform.s_end)];

% Objective Function f*x
f = [ones(1, xform.x_cnt), zeros(1, x_cnt-xform.x_end)];

% linear equality constraint
% set s to zeros
Aeq_set = T_s;
beq_set = zeros(xform.s_cnt, 1);
% infeasible flows
[Aeq_infeasible, beq_infeasible] = flow_infeasible(N, x_cnt, xform, C_r);
Aeq = [Aeq_set; Aeq_infeasible];
beq = [beq_set; beq_infeasible];

% linear inequality constraints
% coverage constraints
Aineq = -cover_matrix(O, N, S_r, x_cnt);
bineq = -repmat(K, N_o, 1);

% bounds
lb = zeros(x_cnt, 1);
ub = [ones(xform.fij_base, 1); Inf(x_cnt-xform.fij_base, 1)];

% integer constraints
nB = xform.x_cnt + xform.s_cnt;           % Number of Binary Variables
nC = xform.fij_cnt + xform.fiB_cnt;       % Number of Continuous Variables
nI = 0;                                   % Number of Integer Variables

% build ctype vector
ctype = [repmat('B', 1, nB), repmat('C', 1, nC), repmat('I', 1, nI)];

% initial guess
x0 = zeros(x_cnt, 1);

% create OPTI Object
%Opt = opti('fun', fun, 'nlmix', nlcon, nlrhs, nle, 'ineq', A, b, ...
%        'bounds', lb, ub, 'eq', Aeq, beq, 'xtype', xtype)

% solve the MINLP problem
%[x,fval,exitflag,info] = solve(Opt,x0)

%options = cplexoptimset('Display', 'on');
[x, fval, exitflag, output] = cplexmilp(f, Aineq, bineq, Aeq, beq, ...
    [], [], [], lb, ub, ctype, x0)

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