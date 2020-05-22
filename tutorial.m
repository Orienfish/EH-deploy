% main.m
clc;
clear;
close all;
warning('off','all');
addpath ('/Applications/CPLEX_Studio1210/cplex/matlab/x86-64_osx');
addpath('./lldistkm');
addpath('./solver');
addpath('./libs');

%% Initialization of the grid map
% pre-process solar and temperature data
fprintf('start pre-processing...\n');
folder = './solardata/';
f_list = dir(append(folder, '*.csv'));
dataT_sav = {append(folder, 'dataT.csv'), append(folder, 'Counts.csv'), ...
    append(folder, 'Centers.csv')}; % all processing data
dataT_sav = string(dataT_sav);
N_bin = 20;        % number of bins
if exist(dataT_sav(1), 'file')
    dataT = readtable(dataT_sav(1));
    Counts = readmatrix(dataT_sav(2));
    Centers = readmatrix(dataT_sav(3));
else
    % if no file exists yet, do pre-process and save it to file
    [dataT, Counts, Centers] = preprocess(f_list, dataT_sav, N_bin);
end

% get data source distribution
bubbleplot_wsize(dataT.lat, dataT.lon, dataT.temp_avg, 'data source temp');
bubbleplot_wsize(dataT.lat, dataT.lon, dataT.dni_avg, 'data source dni');
origin = [min(dataT.lat), min(dataT.lon)];
yScalekm = lldistkm([min(dataT.lat), min(dataT.lon)], [max(dataT.lat), min(dataT.lon)]);
xScalekm = lldistkm([min(dataT.lat), min(dataT.lon)], [min(dataT.lat), max(dataT.lon)]);

% set the size and granularity of the grid space 
xScalem = 1000;    % m
yScalem = 1000;    % m
N_x = 10;
N_y = 10;
N_cnt = N_x * N_y;      % number of grid points
Unit_x = floor(xScalem / (N_x - 1));
Unit_y = floor(yScalem / (N_y - 1));
A = 0.01;               % surface area (m^2) of solar panel
xi = 0.05;              % end-to-end conversion efficiency of solar system

% the factor to transform from grid (m) to lat and lon
% x - longitude, y - latitude
x_transform = (max(dataT.lon) - min(dataT.lon)) / xScalem;
y_transform = (max(dataT.lat) - min(dataT.lat)) / yScalem;

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
        %N(i + j * N_x + 1).Ti = dataT.temp_avg(dataT_idx);
        N(i + j * N_x + 1).Ti = 25 + (40 - 25) * j / N_y;
    end
end
% plot the heatmap of temperature distribution in the grid map
plot_temp(N, N_x, N_y);

%% Initialization of basic parameters
fprintf('Intializing basic parameters...\n');
params.S_r = 120;                       % sensing range in m
params.C_r = 120;                       % communication range in m
params.N_o = 15;                        % number of PoIs
params.K = 1;                           % K-coverage
params.G = 100;                         % generated bytes for each sample
params.B = 2000;                        % bandwidth B/s
params.eta = 0.2;                       % sampling frequency
params.maxf = params.eta*params.G*N_cnt;% maximum flow amount
params.P0 = 0.01;                       % sleep power (W)
params.Es = 0.2*0.2;                    % sensing energy (W*s)
params.Prx = 0.1;                       % reception power (W)
% randomly generate PoIs to monitor
O = repmat([], params.N_o, 2);
for i = 1:params.N_o
    O(i, 1) = unifrnd(0, xScalem);
    O(i, 2) = unifrnd(0, yScalem);
end
% randomly generate the location of the sink
c = [unifrnd(0, xScalem), unifrnd(0, yScalem)];
% get the distance matrix
% dist(i, j) denotes the Euclidean distance between grid i and j
% dist(i, N_cnt+1) denotes the Euclidean distance between i and sink
dist = getDist(N, c);

% call the amb2core function to load the global variables k_1, k_2, k_3
Tcorei = amb2core(25, 3);
% specify the reliability options and targets
rel.SoH = true;
rel.SoHref = 0.8;
rel.T = 5;                              % years
rel.MTTF = true;
rel.MTTFref = 0.75;

% convert the reliability constraints to power constraints
Pi = vertcat(N(:).Ri);                  % power constraints (W)
if rel.SoH == true
    P_soh = Psoh_bound(rel.SoHref, rel.T, vertcat(N(:).Ti));
    Pi = [Pi, P_soh - repmat(params.P0, N_cnt, 1)];
    Pi = max(Pi, 0);
end
if rel.MTTF == true
    P_mttf = Pmttf_bound(rel.MTTFref, vertcat(N(:).Ti));
    Pi = [Pi, P_mttf - repmat(params.P0, N_cnt, 1)];
    Pi = max(Pi, 0);
end
disp(Pi);
Pi = min(Pi, [], 2); % get the column vector of min of each row
for i = 1:N_cnt
    N(i).Pi = Pi(i); % clip the power constraints to grid struct
end

%% Call solvers
% options to run which solver/algorithm
run.cplex = true;
run.tatsh = true;
run.tsh = true;

% Call the CPLEX solver
if run.cplex
    % solve the problem without SoH and MTTF constraints
    rel.SoH = false; rel.MTTF = false;
    sol_wo = solver(N, O, dist, params, rel);
    % plot the solution
    if sol_wo.exitflag == 1
        plot_solution(N, O, c, sol_wo, params.S_r, [xScalem, yScalem]);
        [sol_wo.sohmin, sol_wo.mttfmin, sol_wo.vio] = ...
            rel_check(sol_wo, N, dist, params, rel);
        fprintf('# of nodes of sol_wo: %d\n', sol_wo.fval);
        fprintf('Min SoH: %f Node: %d\n', sol_wo.sohmin(1), sol_wo.sohmin(2));
        fprintf('Min MTTF: %f Node: %d\n', sol_wo.mttfmin(1), sol_wo.mttfmin(2));
        fprintf('Violation of sol_wo: %f\n', sol_wo.vio);
    end
    % solve the problem with SoH and MTTF constraints
    rel.SoH = true; rel.MTTF = true;
    sol_w = solver(N, O, dist, params, rel);
    % plot the solution
    if sol_w.exitflag == 1
        plot_solution(N, O, c, sol_w, params.S_r, [xScalem, yScalem]);
        [sol_w.sohmin, sol_w.mttfmin, sol_w.vio] = ...
            rel_check(sol_w, N, dist, params, rel);
        fprintf('# of nodes of sol_w: %d\n', sol_w.fval);
        fprintf('Min SoH: %f Node: %d\n', sol_w.sohmin(1), sol_w.sohmin(2));
        fprintf('Min MTTF: %f Node: %d\n', sol_w.mttfmin(1), sol_w.mttfmin(2));
        fprintf('Violation of sol_w: %f\n', sol_w.vio);
    end
end

% Call TATSH
tatshparams.w1 = 1;     % weight for placing new node
tatshparams.w2 = 0.05;  % weight for remained power budget
if run.tatsh
    fprintf('calling TATSH...\n');
    sol_tatsh = TATSH(N, O, dist, params, tatshparams);
    % plot the solution
    if sol_tatsh.exitflag == 1
        plot_solution(N, O, c, sol_tatsh, params.S_r, [xScalem, yScalem]);
        [sol_tatsh.sohmin, sol_tatsh.mttfmin, sol_tatsh.vio] = ...
            rel_check(sol_tatsh, N, dist, params, rel);
        fprintf('# of nodes of sol_tatsh: %d\n', sol_tatsh.fval);
        fprintf('Min SoH: %f Node: %d\n', sol_tatsh.sohmin(1), sol_tatsh.sohmin(2));
        fprintf('Min MTTF: %f Node: %d\n', sol_tatsh.mttfmin(1), sol_tatsh.mttfmin(2));
        fprintf('Violation of sol_tatsh: %f\n', sol_tatsh.vio);
    end
end

% Call TSH
tshparams.w1 = 500;     % cost for adding a new node
tshparams.w2 = 800;     % cost for adding per area of solar panel
if run.tatsh
    fprintf('calling TSH...\n');
    sol_tsh = TSH(N, O, dist, params, tatshparams);
    % plot the solution
    if sol_tsh.exitflag == 1
        plot_solution(N, O, c, sol_tsh, params.S_r, [xScalem, yScalem]);
        [sol_tsh.sohmin, sol_tsh.mttfmin, sol_tsh.vio] = ...
            rel_check(sol_tsh, N, dist, params, rel);
        fprintf('# of nodes of sol_tsh: %d\n', sol_tsh.fval);
        fprintf('Min SoH: %f Node: %d\n', sol_tsh.sohmin(1), sol_tsh.sohmin(2));
        fprintf('Min MTTF: %f Node: %d\n', sol_tsh.mttfmin(1), sol_tsh.mttfmin(2));
        fprintf('Violation of sol_tsh: %f\n', sol_tsh.vio);
    end
end

% end of tutorial

%% Appendix functions
% Get distance matrix between grid points
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

% Check the reliability of the solution
% Args:
%   sol: the struct of the solution to be evaluated
%   N: the struct of grid locations
%   dist: the distance matrix between all grid locations and the sink
%   params: necessary basic parameters
%   rel: the struct of reliability options and targets
%
% Return:
%   sohmin: min SoH of all deployed devices
%   mttfmin: min MTTF of all deployed devices
%   vio: percentage of violations among all deployed sites

function [sohmin, mttfmin, vio] = rel_check(sol, N, dist, params, rel)
    % get number of grid locations
    N_cnt = size(N, 1);
    % get the power at all grid locations
    pwr = getPwr(sol, N, dist, params);
    % calculate core temperature
    Tc = zeros(N_cnt, 1);
    for i = 1:N_cnt
        if pwr(i) > 0
            Tc(i) = amb2core(N(i).Ti, pwr(i));
        else
            Tc(i) = 273.15; % in Kelvin, for those undeployed spots
        end
    end
    % calculate minimal SoH of all deployed devices
    SoH = soh(Tc, rel.T);
    sohmin = zeros(1, 2);
    [sohmin(1), sohmin(2)] = min(SoH);
    % calculate minimal MTTF of all deployed devices
    MTTF = mttf(Tc);
    mttfmin = zeros(1, 2);
    [mttfmin(1), mttfmin(2)] = min(MTTF);
    % combine all power bounds from N
    P_bound = vertcat(N(:).Pi);
    % calculate the violation percentage
    vio = sum(pwr > P_bound) / sum(sol.x);    
end

% plot functions
% plot the locations
function bubbleplot_wsize(lat, lon, sizedata, title)
    figure;
    geobubble(lat, lon, sizedata, 'Title', title);
    ax = gca; % get current axes
    ax.FontSize = 16;
    %geobasemap streets-light; % set base map style
end

% plot the heatmap of temperature in the grid space
function plot_temp(N, N_x, N_y)
    temp = flipud(reshape(vertcat(N(:).Ti), [N_x, N_y])');
    figure;
    h = heatmap(1:N_x, 1:N_y, temp, 'Colormap', flipud(autumn), ...
        'CellLabelColor','none', 'FontSize', 16);
end

% plot the solution in the grid space
function plot_solution(N, O, c, sol, S_r, maxlim)
    % intialization
    N_cnt = size(N, 1);  % number of grid locations
    
    % scatter
    color_map = [[0 0.4470 0.7410]; ...       % blue
        [0.8500 0.3250 0.0980]; ...           % orange
        [0.4660 0.6740 0.1880]; ...           % green
        [0.4940 0.1840 0.5560]];              % purple
    nodes = vertcat(N.position);              % grid points
    nodes = vertcat(nodes, O);                % append PoIs
    nodes = vertcat(nodes, c);                % append the sink
    sz_nodes = repmat(40, size(nodes, 1), 1); % const size for nodes
    color_idx = vertcat(sol.x, repmat(2, size(O, 1), 1), 3) + 1;
    color_idx = round(color_idx);
    color_nodes = vertcat(color_map(color_idx, :));
    figure;
    scatter(nodes(:, 1), nodes(:, 2), sz_nodes, color_nodes, 'filled', ...
        'LineWidth', 2);
    hold on;
    
    % plot the coverage circle
    cplot = @(r, x0, y0) plot(x0 + r * cos(linspace(0, 2*pi)), ...
        y0 + r * sin(linspace(0, 2*pi)),'r-', 'LineWidth', 2);
    for i = 1: N_cnt
        if sol.s(i) > 0.5
            cplot(S_r, N(i).position(1), N(i).position(2));
            axis equal; hold on;
        end
    end
    
    % plot the line segment representing flows
    ls = [];             % list of line segments
    for i = 1:N_cnt
        for j = 1:N_cnt
            fij_idx = (i-1) * N_cnt + j;
            if sol.fij(fij_idx) > 0.5 && sol.x(i) > 0.5 && sol.x(j) > 0.5
                ls = [ls; [N(i).position, N(j).position]];
            end
        end
        if sol.fiB(i) > 0.5 && sol.x(i) > 0.5
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
    ax = gca; ax.FontSize = 16;
end