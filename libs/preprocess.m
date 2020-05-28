%% Preprocess the existing data in ./data folder.
% Convert historical readings into average, variance and count.
% A struct with name, location, and statistical parameters are returned.
%
% Args:
%   f_list: list of file names to import data from
%   dataT_sav: the name of files to save extracted dataT, Counts and Centers
%   N_bin: number of bins
%
% Return:
%   dataT: statistical pattern of the data in a table

function [dataT, Counts, Centers] = preprocess(f_list, dataT_sav, N_bin)

% initialize data table dataT
varNames = {'name', 'lat', 'lon', 'obs_cnt', ...
            'temp_avg', 'temp_var', 'temp_cnt', 'temp_max', 'temp_min', ...
            'temp_1', 'temp_5', 'temp_10', 'temp_20', 'temp_50', ...
            'temp_80', 'temp_90', 'temp_95', 'temp_99', ...
            'dni_avg', 'dni_var', 'dni_cnt', 'dni_max', 'dni_min'};
varTypes = {'int32', 'double', 'double', 'int32', ...
            'double', 'double', 'double', 'double', 'double', ...
            'double', 'double', 'double', 'double', 'double', ...
            'double', 'double', 'double', 'double', ...
            'double', 'double', 'int32', 'double', 'double'};
dataT = table('Size', [length(f_list) length(varNames)], ...
              'VariableTypes', varTypes, 'VariableNames', varNames);
% initialize the matrix to store temperature distribution
Counts = zeros(length(f_list), N_bin);
Centers = zeros(length(f_list), N_bin);

figure('Position', [0 0 500 400]);
for f_idx = 1:length(f_list)
    % obtain the list of data files
    f_name = f_list(f_idx).name;
    % fprintf('Processing file %s...\n', f_name);
        
    % split the file name by '()'
    f_split = strsplit(f_name, '_');
    % obtain location id
    dataT.name(f_idx) = str2double(f_split(1));
    % obtain latitude and longtitude
    dataT.lat(f_idx) = str2double(f_split(2));
    dataT.lon(f_idx) = str2double(f_split(3));
    
    % read data
    f_path = append(f_list(f_idx).folder, '/', f_name);
    T = readtable(f_path, 'Delimiter', ',', 'HeaderLines', 2); % jump first two lines
    
    % getthe number of total observations
    dataT.obs_cnt(f_idx) = height(T);
    
    % get temperature in Celsius
    array = T.Temperature(~isnan(T.Temperature)); % filter out nan
    array = array(array >= 0); % filter out outliers
    dataT.temp_avg(f_idx) = mean(array);
    dataT.temp_var(f_idx) = var(array);
    dataT.temp_cnt(f_idx) = length(array);
    dataT.temp_max(f_idx) = max(array);
    dataT.temp_min(f_idx) = min(array);
    
    % fill in the temperature statistical patterns
    dataT.temp_1(f_idx) = quantile(T.Temperature, 0.01);
    dataT.temp_5(f_idx) = quantile(T.Temperature, 0.05);
    dataT.temp_10(f_idx) = quantile(T.Temperature, 0.1);
    dataT.temp_20(f_idx) = quantile(T.Temperature, 0.2);
    dataT.temp_50(f_idx) = quantile(T.Temperature, 0.5);
    dataT.temp_80(f_idx) = quantile(T.Temperature, 0.8);
    dataT.temp_90(f_idx) = quantile(T.Temperature, 0.9);
    dataT.temp_95(f_idx) = quantile(T.Temperature, 0.95);
    dataT.temp_99(f_idx) = quantile(T.Temperature, 0.99);
    
    % get DNI in w/m2
    array = T.DNI(~isnan(T.DNI)); % filter out nan
    array = array(array >= 0); % filter out outliers
    dataT.dni_avg(f_idx) = mean(array);
    dataT.dni_var(f_idx) = var(array);
    dataT.dni_cnt(f_idx) = length(array);
    dataT.dni_max(f_idx) = max(array);
    dataT.dni_min(f_idx) = min(array);
    
    [curCounts, curEdges] = histcounts(T.Temperature, N_bin);
    curCenters = 0.5 * (curEdges(1:N_bin) + curEdges(2:N_bin+1));
    Counts(f_idx, :) = curCounts ./ dataT.temp_cnt(f_idx);
    Centers(f_idx, :) = curCenters;
    % only plot part of the cdfs
    if mod(f_idx, 21) == 0
        h = cdfplot(T.Temperature);
        set(h, 'LineWidth', 1.5);
        hold on;
    end
end
ax = gca; ax.FontSize=16;
xlim([0, 40]);
xlabel('Temperature (°C)'); ylabel('Cumulative Probability');
% save all data to the designated files
writetable(dataT, dataT_sav(1));
writematrix(Counts, dataT_sav(2));
writematrix(Centers, dataT_sav(3));
end