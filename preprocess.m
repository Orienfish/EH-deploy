function dataT = preprocess(f_list, dataT_sav, thres)
% Preprocess the existing data in ./data folder.
% Convert historical readings into average, variance and count.
% A struct with name, location, and statistical parameters are returned.
%
% Args:
%   f_list: list of file names to import data from
%   dataT_sav: the name of file to save extracted dataT
%
% Return:
%   dataT: statistical pattern of the data in a table

% initialize data table dataT
varNames = {'name', 'lat', 'lon', 'obs_cnt', ...
            'temp_avg', 'temp_var', 'temp_cnt', 'temp_max', 'temp_min', ...
            'dni_avg', 'dni_var', 'dni_cnt', 'dni_max', 'dni_min'};
varTypes = {'int32', 'double', 'double', 'int32', ...
            'double', 'double', 'int32', 'double', 'double', ...
            'double', 'double', 'int32', 'double', 'double'};
dataT = table('Size', [length(f_list) length(varNames)], ...
              'VariableTypes', varTypes, 'VariableNames', varNames);

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
    
    % get DNI in w/m2
    array = T.DNI(~isnan(T.DNI)); % filter out nan
    array = array(array >= 0); % filter out outliers
    dataT.dni_avg(f_idx) = mean(array);
    dataT.dni_var(f_idx) = var(array);
    dataT.dni_cnt(f_idx) = length(array);
    dataT.dni_max(f_idx) = max(array);
    dataT.dni_min(f_idx) = min(array);
end
writetable(dataT, dataT_sav);
end