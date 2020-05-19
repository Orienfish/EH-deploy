% main.m
clc;
clear;
close all;
warning('off','all');
addpath ('/Applications/CPLEX_Studio1210/cplex/matlab/x86-64_osx');
addpath('./lldistkm');
addpath('./solver');

%% Initialization of the grid map
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
        %N(i + j * N_x + 1).Ti = dataT.temp_avg(dataT_idx);
        N(i + j * N_x + 1).Ti = 20 + (40 - 20) * j / N_y;
    end
end
% plot the heatmap of temperature distribution in the grid map
plot_temp(N, N_x, N_y);

%% Initialization of basic parameters
fprintf('Intializing basic parameters...\n');
params.S_r = 120;                       % sensing range in m
params.C_r = 120;                       % communication range in m
params.N_o = 10;                        % number of PoIs
params.K = 1;                           % K-coverage
params.G = 100;                         % generated bytes for each sample
params.B = 2000;                        % bandwidth B/s
params.eta = 0.1;                       % sampling frequency
params.maxf = params.eta*params.G*N_cnt;% maximum flow amount
% randomly generate PoIs to monitor
O = repmat([], params.N_o, 2);
for i = 1:params.N_o
    O(i, 1) = unifrnd(0, xScalem_target);
    O(i, 2) = unifrnd(0, yScalem_target);
end
% randomly generate the location of the sink
c = [unifrnd(0, xScalem_target), unifrnd(0, yScalem_target)];
% get the distance matrix
dist = getDist(N, c);

% options to run which solver/algorithm
run.cplex = true;

%% Call the CPLEX solver
if run.cplex == true
    sol = solver(N, O, dist, params);
    % plot the solution
    plot_solution(N, O, c, sol, params.S_r, [xScalem_target, yScalem_target]);
end


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
        if sol.s(i) > 0
            cplot(S_r, N(i).position(1), N(i).position(2));
            axis equal; hold on;
        end
    end
    
    % plot the line segment representing flows
    ls = [];             % list of line segments
    for i = 1:N_cnt
        for j = 1:N_cnt
            fij_idx = (i-1) * N_cnt + j;
            if sol.fij(fij_idx) > 0
                ls = [ls; [N(i).position, N(j).position]];
            end
        end
        if sol.fiB(i) > 0
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