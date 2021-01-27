% main.m
clc;
clear;
close all;
warning('off','all');

%% Experiment 1: small scale simulation
% Set experiment options
% set scale of grid space
exp_opt.xScalem = 500;                         % m
exp_opt.yScalem = 500;                         % m
exp_opt.N_x = 10;
exp_opt.N_y = 10;

% basic parameters
exp_opt.S_r = 100;                       % sensing range in m
exp_opt.C_r = 120;                       % communication range in m
exp_opt.N_o = 20;                        % number of PoIs
exp_opt.K = 2;                           % K-coverage
exp_opt.G = 100;                         % generated bytes for each sample
exp_opt.B = 2000;                        % bandwidth B/s
exp_opt.eta = 0.2;                       % sampling frequency

% specify the reliability options and targets
exp_opt.rel.SoH = true;
exp_opt.rel.SoHref = 0.90;
exp_opt.rel.T = 5;                              % years
exp_opt.rel.MTTF = true;
exp_opt.rel.MTTFref = 0.94;
exp_opt.rel.MTTFsolarref = 1.33;

% options to run which solver/algorithm
run.cplex = true;
run.rdtsh = true;
run.tsh = true;
run.srigh = true;
run.rdsrigh = true;

% which experiment to run
exp.small = true;
exp.large = false;

if exp.small
    % set test rounds
    iter = 20;

    % Call exp functions
    % Subtest 1: various number of targets
    N_o_list = [5, 10, 15, 20, 25];
    res_target_s = [];
    exp_opt.K = 1; % K-coverage
    for i=1:length(N_o_list)
        fprintf('Running small target exp with %d targets\n', N_o_list(i));
        for j=1:iter
            fprintf('iter %d\n', j);
            exp_opt.N_o = N_o_list(i);
            res = exp_func(run, exp_opt);
            while isempty(res)
                res = exp_func(run, exp_opt);
            end
            fill_resT(res, run, './res_target_small.csv');
        end
    end
    exp_opt.N_o = 20; % reset to standard value

    % Subtest 2: various number of sites
    N_x_list = [10, 11, 12, 13, 14];
    res_site_s = [];
    for i=1:length(N_x_list)
        fprintf('Running small site exp with %d sites on x\n', N_x_list(i));
        for j=1:iter
            fprintf('iter %d\n', j);
            exp_opt.N_x = N_x_list(i);
            res = exp_func(run, exp_opt);
            while isempty(res)
                res = exp_func(run, exp_opt);
            end
            fill_resT(res, run, './res_site_small.csv');
        end
    end
    exp_opt.N_x = 10; % reset to standard value
end

%% Experiment 2: large scale simulation
if exp.large
    run.cplex = false;
    % set scale of grid space
    exp_opt.xScalem = 4000;                         % m
    exp_opt.yScalem = 4000;                         % m
    exp_opt.N_x = 80;
    exp_opt.N_y = 80;
    exp_opt.K = 2;                            % K-coverage

    exp_opt.N_o = 100;                        % number of PoIs
    exp_opt.rel.MTTFsolarref = 1;          % solar panel MTTF bound

    % set test rounds
    iter = 20;

    % Call exp functions
    % Subtest 1: various number of targets
    N_o_list = [50, 75, 100, 125, 150];
    res_target_l = [];
    for i=1:length(N_o_list)
        fprintf('Running large target exp with %d targets\n', N_o_list(i));
        for j=1:iter
            fprintf('iter %d\n', j);
            exp_opt.N_o = N_o_list(i);
            res = exp_func(run, exp_opt);
            while isempty(res)
                res = exp_func(run, exp_opt);
            end
            fill_resT(res, run, './res_target_large.csv');
        end
    end
    exp_opt.N_o = 100; % reset to standard value

    % Subtest 2: various number of sites
    N_x_list = [60, 70, 80, 90, 100];
    res_site_l = [];
    for i=1:length(N_x_list)
        fprintf('Running large site exp with %d sites on x\n', N_x_list(i));
        for j=1:iter
            fprintf('iter %d\n', j);
            exp_opt.N_x = N_x_list(i);
            res = exp_func(run, exp_opt);
            while isempty(res)
                res = exp_func(run, exp_opt);
            end
            fill_resT(res, run, './res_site_large.csv');
        end
    end
    exp_opt.N_x = 100; % reset to standard value
end

% Subtest 3: various coverage level
%k_list = [1, 2, 3, 4];
%res_k_l = [];
%for i=1:length(k_list)
%    fprintf('Running large site exp with coverage level %d\n', k_list(i));
%    for j=1:iter
%        fprintf('iter %d\n', j);
%        exp_opt.K = k_list(i);
%        res = exp_func(run, exp_opt);
%        while isempty(res)
%            res = exp_func(run, exp_opt);
%        end
%        res_k_l = [res_k_l; fill_resT(res, run)];
%    end
%end
%exp_opt.K = 2;
%writetable(res_k_l, './res_cover_large.csv');

% Subtest 4: trade-offs between reliability and number of nodes
%run.tsh = false;
%run.srigh = false;
%N_r_list = [0.85, 0.87, 0.89, 0.91, 0.93, 0.95];
%res_trade_l = [];
%for i=1:length(N_r_list)
%    fprintf('Running large target exp rel bound %f\n', N_r_list(i));
%    for j=1:iter
%        fprintf('iter %d\n', j);
%        exp_opt.rel.SoHref = N_r_list(i);
%        exp_opt.rel.MTTFref = N_r_list(i);
%        res = exp_func(run, exp_opt);
%        while isempty(res)
%            res = exp_func(run, exp_opt);
%        end
%        res_trade_l = [res_trade_l; fill_resT(res, run)];
%    end
%end
%exp_opt.rel.SoHref = 0.90;                       % reset to standard value
%exp_opt.rel.MTTFref = 0.90;
%writetable(res_trade_l, './res_trade_large.csv');

%% Appendix functions
% Fill in result table from the result struct and write table to file
% Args:
%   res: the result struct returned by exp_func
%   run: run options
%   filename: the file to write to
%
% Return:
%   resT: one line in the result table
function resT = fill_resT(res, run, filename)
% initialize one line of result table
varNames = {'opt_wo_fval', 'opt_wo_vio', 'opt_wo_sohmin', 'opt_wo_mttfmin', ...
            'opt_wo_time', ...
            'opt_fval', 'opt_vio', 'opt_sohmin', 'opt_mttfmin', ...
            'opt_time', ...
            'rdtsh_fval', 'rdtsh_vio', 'rdtsh_sohmin', 'rdtsh_mttfmin', ...
            'rdtsh_time', ...
            'tsh_fval', 'tsh_vio', 'tsh_sohmin', 'tsh_mttfmin', ...
            'tsh_time', ...
            'rdsrigh_fval', 'rdsrigh_vio', 'rdsrigh_sohmin', 'rdsrigh_mttfmin', ...
            'rdsrigh_time', ...
            'srigh_fval', 'srigh_vio', 'srigh_sohmin', 'srigh_mttfmin', ...
            'srigh_time'};
varTypes = {'double', 'double', 'double', 'double', 'double', ...
            'double', 'double', 'double', 'double', 'double', ...
            'double', 'double', 'double', 'double', 'double', ...
            'double', 'double', 'double', 'double', 'double', ...
            'double', 'double', 'double', 'double', 'double', ...
            'double', 'double', 'double', 'double', 'double'};
resT = table('Size', [1 length(varNames)], ...
    'VariableTypes', varTypes, 'VariableNames', varNames);
% fill in the item one by one if tested
if run.cplex
    % optimization without reliability constraints
    resT.opt_wo_fval(1) = res.sol_wo.fval;
    resT.opt_wo_vio(1) = res.sol_wo.vio;
    resT.opt_wo_sohmin(1) = res.sol_wo.sohmin(1);
    resT.opt_wo_mttfmin(1) = res.sol_wo.mttfmin(1);
    resT.opt_wo_time(1) = res.sol_wo.time;
    % optimization with reliability constraints
    resT.opt_fval(1) = res.sol_w.fval;
    resT.opt_vio(1) = res.sol_w.vio;
    resT.opt_sohmin(1) = res.sol_w.sohmin(1);
    resT.opt_mttfmin(1) = res.sol_w.mttfmin(1);
    resT.opt_time(1) = res.sol_w.time;
end
if run.rdtsh
    resT.rdtsh_fval(1) = res.sol_rdtsh.fval;
    resT.rdtsh_vio(1) = res.sol_rdtsh.vio;
    resT.rdtsh_sohmin(1) = res.sol_rdtsh.sohmin(1);
    resT.rdtsh_mttfmin(1) = res.sol_rdtsh.mttfmin(1);
    resT.rdtsh_time(1) = res.sol_rdtsh.time;
end
if run.tsh
    resT.tsh_fval(1) = res.sol_tsh.fval;
    resT.tsh_vio(1) = res.sol_tsh.vio;
    resT.tsh_sohmin(1) = res.sol_tsh.sohmin(1);
    resT.tsh_mttfmin(1) = res.sol_tsh.mttfmin(1);
    resT.tsh_time(1) = res.sol_tsh.time;
end
if run.rdsrigh
    resT.rdsrigh_fval(1) = res.sol_rdsrigh.fval;
    resT.rdsrigh_vio(1) = res.sol_rdsrigh.vio;
    resT.rdsrigh_sohmin(1) = res.sol_rdsrigh.sohmin(1);
    resT.rdsrigh_mttfmin(1) = res.sol_rdsrigh.mttfmin(1);
    resT.rdsrigh_time(1) = res.sol_rdsrigh.time;
end
if run.srigh
    resT.srigh_fval(1) = res.sol_srigh.fval;
    resT.srigh_vio(1) = res.sol_srigh.vio;
    resT.srigh_sohmin(1) = res.sol_srigh.sohmin(1);
    resT.srigh_mttfmin(1) = res.sol_srigh.mttfmin(1);
    resT.srigh_time(1) = res.sol_srigh.time;
end

% if this is the first line of result, write variables name
if ~exist(filename, 'file')
    writetable(resT, filename);
else % append to existing table
    writetable(resT, filename, 'WriteMode', 'Append', ...
        'WriteVariableNames', false);
end
end