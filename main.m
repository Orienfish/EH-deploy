% main.m
clc;
clear;
close all;
warning('off','all');
addpath('./lldistkm');

%% preparation
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

% get data source D
D_lat = vertcat(dataT.lat(:));
D_lon = vertcat(dataT.lon(:));
D = horzcat(D_lat, D_lon);
bubbleplot_wsize(dataT.lat(:), dataT.lon(:), dataT.temp_avg(:), 'data source');
origin = [min(D_lat), min(D_lon)];
yScalekm = lldistkm([min(D_lat), min(D_lon)], [max(D_lat), min(D_lon)]);
xScalekm = lldistkm([min(D_lat), min(D_lon)], [min(D_lat), max(D_lon)]);

% set the size and granularity of the grid space 
xScalem_target = 2000; % 2000 m
yScalem_target = 2000; % 2000 m
N_x = 10;
N_y = 10;
N_cnt = N_x * N_y;
Unit_x = floor(xScalem_target / N_x);
Unit_y = floor(yScalem_target / N_y);

% the factor to transform from grid (m) to real field (km)
x_transform = xScalekm / xScalem_target;
y_transform = yScalekm / yScalem_target;

% generate the grid candidate set N 
% with their x, y coordinates and temperature and DNI
empty_point.position = [];
empty_point.Ri = [];
empty_point.Ti = [];
N = repmat(empty_point, N_cnt, 1);

for j = 0:N_y-1
    for i = 0:N_x-1
        N(i + j*N_x+1).position = [Unit_x * i, Unit_y * j];
        %x_idx = ;
        %y_idx = ;
        %dataT_idx = ;
        %N(i).Ri = dataT(dataT_idx)....;
        %N(i).Ti = dataT(dataT_idx).temp_avg.
    end
end


%% plot functions
function bubbleplot_wsize(lat, lon, sizedata, title)
    % plot the locations
    figure;
    geobubble(lat, lon, sizedata, 'Title', title);
    ax = gca; % get current axes
    ax.FontSize = 16;
    %geobasemap streets-light; % set base map style
end