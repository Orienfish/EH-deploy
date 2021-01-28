%% Run simulation with specified parameters
%
% Args:
%   run: run which solver/algorithm
%   exp_opt: experiment options

function res = exp_func(run, exp_opt)
warning('off','all');
addpath('/Applications/CPLEX_Studio1210/cplex/matlab/x86-64_osx');
addpath('~/CPLEX_Studio1210/cplex/matlab/x86-64_linux');
addpath('../lldistkm');
addpath('../solver');
addpath('../libs');
addpath('../alg');
addpath('../');
res = [];

%% Initialization of the grid map
% pre-process solar and temperature data
fprintf('start pre-processing...\n');
folder = '../solardata/';
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
%bubbleplot_wsize(dataT.lat, dataT.lon, dataT.temp_avg, 'data source temp');
%bubbleplot_wsize(dataT.lat, dataT.lon, dataT.dni_avg, 'data source dni');
origin = [min(dataT.lat), min(dataT.lon)];
%yScalekm = lldistkm([min(dataT.lat), min(dataT.lon)], [max(dataT.lat), min(dataT.lon)]);
%xScalekm = lldistkm([min(dataT.lat), min(dataT.lon)], [min(dataT.lat), max(dataT.lon)]);

% set the size and granularity of the grid space
xScalem = exp_opt.xScalem;    % m
yScalem = exp_opt.yScalem;    % m
N_x = exp_opt.N_x;
N_y = exp_opt.N_y;
N_cnt = N_x * N_y;      % number of grid points
Unit_x = floor(xScalem / (N_x - 1));
Unit_y = floor(yScalem / (N_y - 1));
A = 0.01;               % surface area (m^2) of solar panel

% the factor to transform from grid (m) to lat and lon
% x - longitude, y - latitude
x_transform = (max(dataT.lon) - min(dataT.lon)) / xScalem;
y_transform = (max(dataT.lat) - min(dataT.lat)) / yScalem;

% generate the grid candidate set N
% with their x, y coordinates and temperature and DNI
empty_point.position = [];
empty_point.Ri = [];    % recharging power after efficiency conversion
empty_point.Ti = [];    % temperature
empty_point.xi = [];    % end-to-end conversion efficiency
empty_point.MTTFi = []; % MTTF of solar panel in ratio
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
        % memorize the original idx for final exportation
        N(i + j * N_x + 1).dataT_idx = dataT_idx;
        % assign the corresponding temperature distribution in Celsius
        N(i + j * N_x + 1).Ti = dataT.temp_avg(dataT_idx) + 4.0;
        N(i + j * N_x + 1).Tcen = Centers(dataT_idx, :) + 4.0;
        N(i + j * N_x + 1).Tcnt = Counts(dataT_idx, :);
        % assign conversion efficiency according to average temperature
        N(i + j * N_x + 1).xi = eff(N(i + j * N_x + 1).Ti, A, ...
            dataT.dni_avg(dataT_idx));
        % assign the corresponding recharging power in W
        % conversion efficiency, solar panel area, w/m2 radiation
        N(i + j * N_x + 1).Ri = N(i + j * N_x + 1).xi * A * ...
            dataT.dni_avg(dataT_idx);
        % compute the expectation of MTTF of solar panel in ratio
        N(i + j * N_x + 1).MTTFi = mttf_solar(N(i + j * N_x + 1));
    end
end
% plot the heatmap of temperature distribution in the grid map
%plot_temp(N, N_x, N_y);

%% Initialization of basic parameters
fprintf('Initializing basic parameters...\n');
params.S_r = exp_opt.S_r;               % sensing range in m
params.C_r = exp_opt.C_r;               % communication range in m
params.N_o = exp_opt.N_o;               % number of PoIs
params.K = exp_opt.K;                   % K-coverage
params.G = 100;                         % generated bytes for each sample
params.B = 2000;                        % bandwidth B/s
params.eta = exp_opt.eta;               % sampling frequency
params.maxf = params.eta*params.G*N_cnt;% maximum flow amount
params.P0 = 0.01;                       % sleep power (W)
params.Es = 0.2*0.2;                    % sensing energy (W*s)
params.Prx = 0.1;                       % reception power (W)

% specify the reliability options and targets
%rel.SoH = true;
%rel.SoHref = 0.90;
%rel.T = 5;                              % years
%rel.MTTF = true;
%rel.MTTFref = 0.90;
%rel.MTTFsolarref = 1.33;
rel = exp_opt.rel;

% eliminate the positions that violate the solar panel reliability bound
MTTFi = vertcat(N(:).MTTFi);
N = N(MTTFi > rel.MTTFsolarref);
N_cnt = size(N, 1);

% convert the reliability constraints to power constraints
Pi = vertcat(N(:).Ri) ;                 % power constraints (W)
if rel.SoH == true
    % use piece-wise approximation of ambient temperature over time
    P_soh = Psoh_bound(rel, N);
    % use average ambient temperation at one location
    %P_soh_Tavg = Psoh_bound_Tavg(rel.SoHref, rel.T, vertcat(N(:).Ti));
    Pi = [Pi, P_soh];
end
if rel.MTTF == true
    % use piece-wise approximation of ambient temperature over time
    P_mttf = Pmttf_bound(rel, N);
    % use average ambient temperation at one location
    %P_mttf_Tavg = Pmttf_bound_Tavg(rel.MTTFref, vertcat(N(:).Ti));
    Pi = [Pi, P_mttf];
end
%disp(Pi);
Pi = min(Pi, [], 2); % get the column vector of min of each row
for idx = 1:N_cnt
    N(idx).Pi = Pi(idx); % clip the power constraints to grid struct
end

% randomly generate PoIs to monitor
O = repmat([], params.N_o, 2);
for i = 1:params.N_o
    flag = params.K;   % how many nodes can cover this target
    while flag > 0
        flag = params.K;
        O(i, 1) = unifrnd(0, xScalem);
        O(i, 2) = unifrnd(0, yScalem);
        for j = 1:N_cnt
            if norm(O(i, :) - N(j).position) < params.S_r
                flag = flag - 1;
            end
            if flag <= 0
                break; % early termination
            end
        end
    end
end

% randomly generate the location of the sink
c = [unifrnd(0, xScalem), unifrnd(0, yScalem)];
% get the distance matrix
% dist(i, j) denotes the Euclidean distance between grid i and j
% dist(i, N_cnt+1) denotes the Euclidean distance between i and sink
dist = getDist(N, c);

%% Call solvers
% options to run which solver/algorithm
%run.cplex = true;
%run.rdtsh = true;
%run.tsh = true;
%run.srigh = true;
%run.rdsrigh = true;


% Call the CPLEX solver
if run.cplex
    % solve the problem with SoH and MTTF constraints
    rel.SoH = true; rel.MTTF = true;
    tic
    sol_w = solver(N, O, dist, params, rel);
    sol_w.time = toc;
    % plot the solution
    if sol_w.exitflag == 1
        %plot_solution(N, O, c, sol_w, params.S_r, ...
        %    [xScalem, yScalem], 'CPLEX w/ Rel');
        sol_w = rel_check(sol_w, N, dist, params, rel);
        log('OPT', sol_w);
        %export_solution(N, c, sol_w, dist, dataT, 'OPT');
    else
        fprintf('No feasible solution for OPT_w!\n');
        return;
    end
    res.sol_w = sol_w;
    % solve the problem without SoH and MTTF constraints
    rel.SoH = false; rel.MTTF = false;
    tic
    sol_wo = solver(N, O, dist, params, rel);
    sol_wo.time = toc;
    % plot the solution
    if sol_wo.exitflag == 1
        %plot_solution(N, O, c, sol_wo, params.S_r, ...
        %    [xScalem, yScalem], 'CPLEX w/o Rel');
        sol_wo = rel_check(sol_wo, N, dist, params, rel);
        log('OPT_wo', sol_wo);
        %export_solution(N, c, sol_wo, dist, dataT, 'OPT_wo');
    else
        fprintf('No feasible solution for OPT_wo!\n');
        return;
    end
    res.sol_wo = sol_wo;
end

% Call RDTSH
if run.rdtsh
    fprintf('calling RDTSH...\n');
    rdtshparams.w1 = 1;  % weight for placing new node
    rdtshparams.w2 = 1.5;  % weight for remained power budget
    tic
    sol_rdtsh = RDTSH(N, O, dist, params, rdtshparams);
    sol_rdtsh.time = toc;
    % plot the solution
    if sol_rdtsh.exitflag == 1
        %plot_solution(N, O, c, sol_rdtsh, params.S_r, ...
        %    [xScalem, yScalem], 'RDTSH');
        sol_rdtsh = rel_check(sol_rdtsh, N, dist, params, rel);
        log('RDTSH', sol_rdtsh);
        %export_solution(N, c, sol_rdtsh, dist, dataT, 'RDTSH');
    else
        fprintf('No feasible solution for RDTSH!\n');
        return;
    end
    res.sol_rdtsh = sol_rdtsh;  
end

% Call TSH
if run.tsh
    fprintf('calling TSH...\n');
    tshparams.w1 = 1;     % cost for adding a new node
    tshparams.w2 = 1;     % cost for adding per area of solar panel
    tic
    sol_tsh = TSH(N, O, dist, params, tshparams);
    sol_tsh.time = toc;
    % plot the solution
    if sol_tsh.exitflag == 1
        %plot_solution(N, O, c, sol_tsh, params.S_r, ...
        %    [xScalem, yScalem], 'TSH');
        sol_tsh = rel_check(sol_tsh, N, dist, params, rel);
        log('TSH', sol_tsh);
        %export_solution(N, c, sol_tsh, dist, dataT, 'TSH');
    else
        fprintf('No feasible solution for TSH!\n');
        return;
    end
    res.sol_tsh = sol_tsh;
end

% Call RDSRIGH
if run.rdsrigh
    fprintf('calling RDSRIGH...\n');
    rdsrighparams.w1 = 1;      % cost for adding a new node
    rdsrighparams.w2 = 12;     % cost for adding per area of solar panel
    tic
    sol_rdsrigh = RDSRIGH(N, O, dist, params, rdsrighparams);
    sol_rdsrigh.time = toc;
    if sol_rdsrigh.exitflag == 1
        %plot_solution(N, O, c, sol_rdsrigh, params.S_r, ...
        %    [xScalem, yScalem], 'RDSRIGH');
        sol_rdsrigh = rel_check(sol_rdsrigh, N, dist, params, rel);
        log('RDSRIGH', sol_rdsrigh);
        %export_solution(N, c, sol_rdsrigh, dist, dataT, 'RDSRIGH');
    else
        fprintf('No feasible solution for RDSRIGH!\n');
        return;
    end
    res.sol_rdsrigh = sol_rdsrigh;
end

% Call SRIGH
if run.srigh
    fprintf('calling SRIGH...\n');
    srighparams.w1 = 1;      % cost for adding a new node
    srighparams.w2 = 1;      % cost for adding per area of solar panel
    tic
    sol_srigh = SRIGH(N, O, dist, params, srighparams);
    sol_srigh.time = toc;
    if sol_srigh.exitflag == 1
        %plot_solution(N, O, c, sol_srigh, params.S_r, ...
        %    [xScalem, yScalem], 'SRIGH');
        sol_srigh = rel_check(sol_srigh, N, dist, params, rel);
        log('SRIGH', sol_srigh);
        %export_solution(N, c, sol_srigh, dist, dataT, 'SRIGH');
    else
        fprintf('No feasible solution for SRIGH!\n');
    end
    res.sol_srigh = sol_srigh;
end

% end of exp func
end

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
%   sol: updated sol struct with SoH and MTTF and violation percentage
function [sol] = rel_check(sol, N, dist, params, rel)
    N_cnt = size(N, 1);         % get number of grid locations
    N_bin = length(N(1).Tcen);  % number of temperature bins
    % get the power at all grid locations
    pwr = getPwr(sol, N, dist, params);
    % calculate core temperature
    Tc = zeros(N_cnt, 1);
    SoH = zeros(N_cnt, 1);
    MTTF = zeros(N_cnt, 1);
    for i = 1:N_cnt
        if pwr(i) > 0
            for j = 1:N_bin
                Tc(i) = amb2core(N(i).Tcen(j), pwr(i));
                SoH(i) = SoH(i) + N(i).Tcnt(j) * soh(Tc(i), rel.T);
                MTTF(i) = MTTF(i) + N(i).Tcnt(j) * mttf(Tc(i));
            end
        else
            Tc(i) = 273.15; % in Kelvin, for those undeployed spots
            SoH(i) = soh(Tc(i), rel.T);
            MTTF(i) = mttf(Tc(i));
        end
    end
    % calculate minimal SoH of all deployed devices
    sol.SoH = SoH;
    sol.sohmin = zeros(1, 2);
    [sol.sohmin(1), sol.sohmin(2)] = min(SoH);
    
    % calculate minimal MTTF of all deployed devices
    sol.MTTF = MTTF;
    sol.mttfmin = zeros(1, 2);
    [sol.mttfmin(1), sol.mttfmin(2)] = min(MTTF);
    
    % combine all power bounds from N
    %P_bound = vertcat(N(:).Pi);
    % calculate the violation percentage
    sol.vio = (sum(SoH < rel.SoHref) + sum(MTTF < rel.MTTFref)) / ...
        (2 * sum(sol.x));
end

% logging function
function log(name, sol)
    fprintf('# of nodes of %s: %d\n', name, sol.fval);
    fprintf('Min SoH: %f Node: %d\n', sol.sohmin(1), sol.sohmin(2));
    fprintf('Min MTTF: %f Node: %d\n', sol.mttfmin(1), sol.mttfmin(2));
    fprintf('Violation of %s: %f\n', name, sol.vio);
    fprintf('Execution time of %s: %f\n', name, sol.time);
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
 %(P_inc(i, j)/max(1e-10, (N(i).Pi - params.P0)))
% plot the heatmap of temperature in the grid space
function plot_temp(N, N_x, N_y)
    temp = flipud(reshape(vertcat(N(:).Ti), [N_x, N_y])');
    figure;
    h = heatmap(1:N_x, 1:N_y, temp, 'Colormap', flipud(autumn), ...
        'CellLabelColor','none', 'FontSize', 16);
end

% plot the solution in the grid space
function plot_solution(N, O, c, sol, S_r, maxlim, method)
    % intialization
    N_cnt = size(N, 1);  % number of grid locations

    % scatter
    color_map = [[0 0.4470 0.7410]; ...       % blue
        [0.3010 0.7450 0.9330]; ...           % ryan
        [0.8500 0.3250 0.0980]; ...           % orange
        [0.4660 0.6740 0.1880]; ...           % green
        [0.6350 0.0780 0.1840]];              % red
    nodes = vertcat(N.position);              % grid points
    nodes = vertcat(nodes, O);                % append PoIs
    nodes = vertcat(nodes, c);                % append the sink
    sz_nodes = repmat(80, size(nodes, 1), 1); % const size for nodes
    color_idx = vertcat(sol.x+sol.s, repmat(3, size(O, 1), 1), 4) + 1;
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
    xlabel('x (m)'); ylabel('y (m)');
    ax = gca; ax.FontSize = 16;
    title(method);
end

% export the solution to text file
function export_solution(N, c, sol, dist, dataT, method)
    N_cnt = size(N, 1);         % get number of grid locations
    
    % create one result folder if it doesn't exist
    if ~exist('res', 'dir')
       mkdir('res')
    end
    
    % export the deployed sensors and corresponding temperature traces
    % of every 4 hours at the real-world location
    filename = sprintf('res/sr_%d_%s.txt', floor(sol.fval), method);
    fileID = fopen(filename, 'w');
    filetemp = sprintf('res/temp_%d_%s.txt', floor(sol.fval), method);
    filetempID = fopen(filetemp, 'w');
    for i=1:N_cnt
        if sol.x(i) > 0.5
            fprintf(fileID, '%.2f %.2f %d\n', N(i).position(1), ...
                N(i).position(2), sol.s(i));
            % search for the corresponding original trace file
            dataT_idx = N(i).dataT_idx;
            pos_str = sprintf('%.2f_%.2f', dataT.lat(dataT_idx), ...
                dataT.lon(dataT_idx));
            folder = '../solardata/';
            f_list = dir(append(folder, sprintf('*%s*.csv', pos_str)));
            if isempty(f_list)
                fprintf(['Solution Export Error! No trace file!\n', ...
                    'Please make sure you download the orignal dataset ', ...
                    'following the README instructions!\n']);
                return;
            end
            f_name = f_list(1).name;
            T = readtable(append(folder, f_name), 'Delimiter', ',', ...
                'HeaderLines', 2); % jump first two lines
            array = T.Temperature(~isnan(T.Temperature)); % filter out nan
            for t_idx=1:floor(size(array, 1)/8)
                fprintf(filetempID, '%.2f ', mean(array(8*(t_idx-1)+1:8*t_idx)));
            end
            fprintf(filetempID, '\n');
        end
    end
    fclose(fileID);
    fclose(filetempID);
    
    % export the deployed sink
    filename = sprintf('res/gw_%d_%s.txt', floor(sol.fval), method);
    fileID = fopen(filename, 'w');
    fprintf(fileID, '%f %f\n', c(1), c(2));
    fclose(fileID);
    
    % export the flow matrix
    filename = sprintf('res/fl_%d_%s.txt', floor(sol.fval), method);
    fileID = fopen(filename, 'w');
    for i=1:N_cnt
        % export flow to other relay nodes
        if sol.x(i) > 0.5
            for j=1:N_cnt
                if sol.x(j) > 0.5
                    fij_idx = (i-1) * N_cnt + j;
                    fprintf(fileID, '%.2f ', sol.fij(fij_idx));
                end
            end
            % export flow to the sink
            fprintf(fileID, '%.2f\n', sol.fiB(i)); 
        end
    end
    fclose(fileID);
    
    % export the transmission power matrix
    filename = sprintf('res/ptx_%d_%s.txt', floor(sol.fval), method);
    fileID = fopen(filename, 'w');
    for i=1:N_cnt
        % only consider placed nodes
        if sol.x(i) > 0.5
            fij_array = sol.fij((i-1)*N_cnt+1 : i*N_cnt);
            % add the last flag for node-sink connection
            flag = vertcat(fij_array > 0.5, sol.fiB(i) > 0.5);
            dist_array = dist(i, logical(flag));
            % compute the max transmission distance
            % export the max transmission power of each placed node
            [max_dist, max_idx] = max(dist_array);
            fprintf(fileID, '%.2f\n', getPtx(max_dist));
        end
    end
    fclose(fileID);
end
