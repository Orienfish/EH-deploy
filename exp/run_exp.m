% main.m
clc;
clear;
close all;
warning('off','all');

%% Experiment 1: small scale simulation
% Set experiment options
% set scale of grid space
exp_opt.xScalem = 1000;                         % m
exp_opt.yScalem = 1000;                         % m
exp_opt.N_x = 10;
exp_opt.N_y = 10;

% basic parameters
exp_opt.S_r = 120;                       % sensing range in m
exp_opt.C_r = 200;                       % communication range in m
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
exp_opt.rel.MTTFref = 0.90;

% options to run which solver/algorithm
run.cplex = true;
run.rdtsh = true;
run.tsh = true;
run.srigh = true;

% set test rounds
iter = 40;

% Call exp functions
% Subtest 1: various number of targets
%N_o_list = [5, 10, 15, 20, 25, 30];
%res_target_s = [];
%for i=1:length(N_o_list)
%    fprintf('Running small target exp with %d targets\n', N_o_list(i));
%    for j=1:iter
%        fprintf('iter %d\n', j);
%        exp_opt.N_o = N_o_list(i);
%        res = exp_func(run, exp_opt);
%        while isempty(res)
%            res = exp_func(run, exp_opt);
%        end
%        res_target_s = [res_target_s; fill_resT(res, run)];
%    end
%end
%exp_opt.N_o = 20;                       % reset to standard value
%writetable(res_target_s, './res_target_small.csv');

% Subtest 2: various number of sites
%N_x_list = [9, 10, 11, 12, 13];
%res_site_s = [];
%for i=1:length(N_x_list)
%    fprintf('Running small site exp with %d sites on x\n', N_x_list(i));
%    for j=1:iter
%        fprintf('iter %d\n', j);
%        exp_opt.N_x = N_x_list(i);
%        res = exp_func(run, exp_opt);
%        while isempty(res)
%            res = exp_func(run, exp_opt);
%        end
%        res_site_s = [res_site_s; fill_resT(res, run)];
%    end
%end
%exp_opt.N_x = 10;
%writetable(res_site_s, './res_site_small.csv');

%% Experiment 2: large scale simulation
run.cplex = false;
% set scale of grid space
exp_opt.xScalem = 10000;                         % m
exp_opt.yScalem = 10000;                         % m
exp_opt.N_x = 100;
exp_opt.N_y = 100;

exp_opt.N_o = 100;                        % number of PoIs

% set test rounds
iter = 40;

% Call exp functions
% Subtest 1: various number of targets
N_o_list = [50, 75, 100, 125, 150, 175];
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
        res_target_l = [res_target_l; fill_resT(res, run)];
    end
end
exp_opt.N_o = 100;                       % reset to standard value
writetable(res_target_l, './res_target_large.csv');

% Subtest 2: various number of sites
N_x_list = [70, 80, 90, 100, 110, 120];
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
        res_site_l = [res_site_l; fill_resT(res, run)];
    end
end
exp_opt.N_x = 100;
writetable(res_site_l, './res_site_large.csv');

%% Appendix functions
% Fill in result table from the result struct
% Args:
%   res: the result struct returned by exp_func
%   run: run options
%
% Return:
%   resT: one line in the result table
function resT = fill_resT(res, run)
% initialize one line of result table
varNames = {'opt_wo_fval', 'opt_wo_vio', 'opt_wo_sohmin', 'opt_wo_mttfmin', ...
            'opt_wo_time', ...
            'opt_fval', 'opt_vio', 'opt_sohmin', 'opt_mttfmin', ...
            'opt_time', ...
            'rdtsh_fval', 'rdtsh_vio', 'rdtsh_sohmin', 'rdtsh_mttfmin', ...
            'rdtsh_time', ...
            'tsh_fval', 'tsh_vio', 'tsh_sohmin', 'tsh_mttfmin', ...
            'tsh_time', ...
            'srigh_fval', 'srigh_vio', 'srigh_sohmin', 'srigh_mttfmin', ...
            'srigh_time'};
varTypes = {'double', 'double', 'double', 'double', 'double', ...
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
if run.srigh
    resT.srigh_fval(1) = res.sol_srigh.fval;
    resT.srigh_vio(1) = res.sol_srigh.vio;
    resT.srigh_sohmin(1) = res.sol_srigh.sohmin(1);
    resT.srigh_mttfmin(1) = res.sol_srigh.mttfmin(1);
    resT.srigh_time(1) = res.sol_srigh.time;
end
end