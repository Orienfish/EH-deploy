% main.m
clc;
clear;
close all;
warning('off','all');
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
xScalem_target = 1000;      % 2000 m
yScalem_target = 1000;      % 2000 m
N_x = 11;
N_y = 11;
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
S_r = 100;          % sensing range in m
C_r = 120;          % communication range in m
N_o = 15;           % number of PoIs
% randomly generate PoIs to monitor
O = repmat([], N_o, 1);
for i = 1:N_o
    O(i, 1) = unifrnd(0, xScalem_target);
    O(i, 2) = unifrnd(0, yScalem_target);
end
% coverage constraint
A = cover_matrix(O, N, S_r);

%% call the solver
% objective
fun = @(x) ones(1, N_cnt)*x;    % Objective Function f(x)

% constraints
A = -A;                         % Linear Inequality Constraints (Ax <= b)
b = -ones(N_o, 1);

% bounds
lb = zeros(N_cnt, 1);
ub = ones(N_cnt, 1);

% integer constraints
nC = 0;                        % Number of Continuous Variables
nI = 0;                        % Number of Integer Variables
nB = N_cnt;                        % Number of Binary Variables

% Build xtype vector
xtype = [repmat('C', 1, nC), repmat('I', 1, nI), repmat('B', 1, nB)];

% initial guess
x0 = zeros(N_cnt, 1);

% Create OPTI Object
%Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'ineq',A,b,'bounds',lb,ub,...
%           'xtype',xtype)
Opt = opti('fun', fun, 'ineq', A, b, 'bounds', lb, ub, 'xtype', xtype)

% Solve the MINLP problem
[x,fval,exitflag,info] = solve(Opt,x0)

%% plot functions
function bubbleplot_wsize(lat, lon, sizedata, title)
    % plot the locations
    figure;
    geobubble(lat, lon, sizedata, 'Title', title);
    ax = gca; % get current axes
    ax.FontSize = 16;
    %geobasemap streets-light; % set base map style
end