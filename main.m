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
bubbleplot_wsize(dataT.lat, dataT.lon, dataT.temp_avg, 'data source temp');
bubbleplot_wsize(dataT.lat, dataT.lon, dataT.dni_avg, 'data source dni');
origin = [min(dataT.lat), min(dataT.lon)];
yScalekm = lldistkm([min(dataT.lat), min(dataT.lon)], [max(dataT.lat), min(dataT.lon)]);
xScalekm = lldistkm([min(dataT.lat), min(dataT.lon)], [min(dataT.lat), max(dataT.lon)]);

% set the size and granularity of the grid space 
xScalem_target = 1000;    % m
yScalem_target = 1000;    % m
N_x = 10;
N_y = 10;
N_cnt = N_x * N_y;      % number of grid points
Unit_x = floor(xScalem_target / (N_x - 1));
Unit_y = floor(yScalem_target / (N_y - 1));
A = 0.01;               % surface area (m^2) of solar panel
xi = 0.05;              % end-to-end conversion efficiency of solar system

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

% project the online dataset to our grid space
for j = 0:N_y-1
    for i = 0:N_x-1
        N(i + j * N_x + 1).position = [Unit_x * i, Unit_y * j];
        x_lon = Unit_x * i * x_transform + origin(2);
        y_lat = Unit_y * j * y_transform + origin(1);
        % find the closest location in dataT
        [minValue, dataT_idx] = min(abs(dataT.lat - y_lat) + ...
            abs(dataT.lon - x_lon));
        % assign the corresponding recharging power in W
        % conversion efficiency, solar panel area, w/m2 radiation
        N(i + j * N_x + 1).Ri = xi * A * dataT.dni_avg(dataT_idx);
        % assign the corresponding temperature in Celsius
        N(i + j * N_x + 1).Ti = dataT.temp_avg(dataT_idx);
    end
end

%% prepare for constraints
fprintf('Preparing for constraints...\n');
S_r = 120;           % sensing range in m
C_r = 120;           % communication range in m
N_o = 5;            % number of PoIs
K = 1;              % K-coverage
G = 100;            % generated bytes for each sample
B = 2000;           % bandwidth B/s
eta = 0.1;          % sampling frequency
maxf = eta*G*N_cnt; % maximum flow amount
% randomly generate PoIs to monitor
O = repmat([], N_o, 2);
for i = 1:N_o
    O(i, 1) = unifrnd(0, xScalem_target);
    O(i, 2) = unifrnd(0, yScalem_target);
end
% location of the sink
c = [xScalem_target/2, yScalem_target/2];
% get the distance matrix
dist = getDist(N, c);

% the binary vectors to locate each type of variable
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
[Afc, bfc] = flow_conserve(G, eta, v_cnt, xform, dist, C_r);
% comple connectivity constraints
v_fiB = [zeros(1, xform.fiB_base), ones(1, xform.fiB_cnt)];
v_s = [zeros(1, xform.x_cnt), ones(1, xform.s_cnt), ...
    zeros(1, v_cnt - xform.s_end)];
Acc = v_fiB - eta * G * v_s; % 1*1
bcc = 0;
% stack all equality constraints
Aeq = [Afc; Acc];
beq = [bfc; bcc];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear inequality constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coverage constraints
[Ac, bc] = coverage(O, N, S_r, v_cnt, K);
% xi >= si
% get the transformation matrix to extract x and s
T_x = [eye(xform.x_cnt), zeros(xform.x_cnt, v_cnt-xform.x_end)];
T_s = [zeros(xform.s_cnt, xform.x_cnt), eye(xform.s_cnt), ...
    zeros(xform.s_cnt, v_cnt-xform.s_end)];
Axs = T_s - T_x;
bxs = zeros(N_cnt, 1);
% sigma_j fij <= maxf * xi
[Abound, bbound] = flow_bound(v_cnt, xform, dist, C_r, maxf);
% power inequality
[Apwr, bpwr] = power_eno(eta, B, v_cnt, xform, N, dist, C_r);
% stack all inequality constraints
Aineq = [Ac; Axs; Abound; Apwr];
bineq = [bc; bxs; bbound; bpwr];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lb = zeros(v_cnt, 1);
ub = [ones(xform.fij_base, 1); repmat(maxf, v_cnt-xform.fij_base, 1)];

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
options.Display = 'On';
[x, fval, exitflag, output] = cplexmilp(f, Aineq, bineq, Aeq, beq, ...
    [], [], [], lb, ub, ctype, x0)

% print the solution
%T_x * x
%T_s * x
% extract flow matrix fij (N_cnt * N_cnt) from x
%T_fij = [zeros(xform.fij_cnt, xform.fij_base), eye(xform.fij_cnt), ...
%    zeros(xform.fij_cnt, v_cnt - xform.fij_end)];
%fij = T_fij * x;
%fij = reshape(fij, N_cnt, N_cnt)'; fij % (i, j) is flow from i to j 
% extract flow to sink fiB (N_cnt * 1) from x
T_fiB = [zeros(xform.fiB_cnt, xform.fiB_base), eye(xform.fiB_cnt)];
fiB = T_fiB * x; fiB                   % (i) is flow from i to B

% plot the solution
plot_solution(N, O, c, x, xform, S_r, [xScalem_target, yScalem_target]);

%% Get distance matrix between grid points
% Args:
%   N: the struct of grid locations
%   c: the location of the sink
%
% Return:
%   dist: the N_cnt * (N_cnt+1) distance matrix.
%         dist(i, j) denotes the Euclidean distance between grid i and j
%         dist(i, N_cnt+1) denotes the Euclidean distance between i and
%         sink
function dist = getDist(N, c)
N_cnt = size(N, 1);
dist = zeros(N_cnt, N_cnt+1);
for i = 1:N_cnt
    for j = i+1:N_cnt
        d = norm(N(i).position - N(j).position);
        dist(i, j) = d;
        dist(j, i) = d;
    end
    dist(i, N_cnt+1) = norm(N(i).position - c);
end
end


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
function plot_solution(N, O, c, sol, xform, S_r, maxlim)
    % intialization
    N_cnt = size(N, 1);  % number of grid locations
    
    % extract different portion of the answer from solution vector
    x = sol(xform.x_base+1: xform.x_end);
    s = sol(xform.s_base+1: xform.s_end);
    fij = sol(xform.fij_base+1: xform.fij_end);
    fiB = sol(xform.fiB_base+1: xform.fiB_end);
    
    % scatter
    color_map = [[0 0.4470 0.7410]; ...       % blue
        [0.8500 0.3250 0.0980]; ...           % orange
        [0.4660 0.6740 0.1880]; ...           % green
        [0.4940 0.1840 0.5560]];              % purple
    nodes = vertcat(N.position);              % grid points
    nodes = vertcat(nodes, O);                % append PoIs
    nodes = vertcat(nodes, c);                % append the sink
    sz_nodes = repmat(40, size(nodes, 1), 1); % const size for nodes
    color_idx = vertcat(x, repmat(2, size(O, 1), 1), 3) + 1;
    color_nodes = vertcat(color_map(color_idx, :));
    figure;
    scatter(nodes(:, 1), nodes(:, 2), sz_nodes, color_nodes, 'filled', ...
        'LineWidth', 2);
    hold on;
    
    % plot the coverage circle
    cplot = @(r, x0, y0) plot(x0 + r * cos(linspace(0, 2*pi)), ...
        y0 + r * sin(linspace(0, 2*pi)),'r-', 'LineWidth', 2);
    for i = 1: N_cnt
        if s(i) > 0
            cplot(S_r, N(i).position(1), N(i).position(2));
            axis equal; hold on;
        end
    end
    
    % plot the line segment representing flows
    ls = [];             % list of line segments
    for i = 1:N_cnt
        for j = 1:N_cnt
            fij_idx = (i-1) * N_cnt + j;
            if fij(fij_idx) > 0
                ls = [ls; [N(i).position, N(j).position]];
            end
        end
        if fiB(i) > 0
            ls = [ls; [N(i).position, c]];
        end
    end
    % plot all line segments
    for i = 1:size(ls, 1)
        plot([ls(i, 1), ls(i, 3)], [ls(i, 2), ls(i, 4)], 'b-', ...
            'LineWidth', 2);
        hold on;
    end
    xlim([0, maxlim(1)]); ylim([0, maxlim(2)]);
end