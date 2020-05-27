%% plot functions for the paper
clc;
clear;
close all;
warning('off','all');
iter = 40;

% small-scale simulation with various targets
res_target_s = readtable('./res_target_small.csv');
N_o_list = [5, 10, 15, 20, 25, 30];
% initialize struct
OPT = [];
OPTnr = [];
RDTSH = [];
TSH = [];
SRIGH = [];
% fill in exp results
for i=1:length(N_o_list)
    [OPTnew, OPTnrnew, RDTSHnew, TSHnew, SRIGHnew] = ...
        fill_res(res_target_s, i, iter);
    OPT = [OPT; OPTnew];
    OPTnr = [OPTnr; OPTnrnew];
    RDTSH = [RDTSH; RDTSHnew];
    TSH = [TSH; TSHnew];
    SRIGH = [SRIGH; SRIGHnew];
end
% plot
figure;
subplot(1, 2, 1);
plot(N_o_list, vertcat(OPT(:).fval), '-*', 'LineWidth', 2, 'DisplayName', 'OPT'); hold on;
plot(N_o_list, vertcat(OPTnr(:).fval),'-o', 'LineWidth', 2, 'DisplayName', 'OPTnr'); hold on;
plot(N_o_list, vertcat(RDTSH(:).fval), '-s', 'LineWidth', 2, 'DisplayName', 'RDTSH'); hold on;
plot(N_o_list, vertcat(TSH(:).fval), '-d', 'LineWidth', 2, 'DisplayName', 'TSH'); hold on;
plot(N_o_list, vertcat(SRIGH(:).fval), '-^', 'LineWidth', 2, 'DisplayName', 'SRIGH'); hold on;
legend(); xlim([5, 30]);
xlabel('Number of PoIs'); ylabel('Number of Nodes');
ax = gca; ax.FontSize = 16;


% fill in the struct with mean value in the 40 iterations
function [OPT, OPTnr, RDTSH, TSH, SRIGH] = fill_res(table, i, iter)
    st_idx = (i-1)*iter+1;
    ed_idx = i*iter;
    OPT.fval = mean(table.opt_fval(st_idx:ed_idx));
    OPT.vio = mean(table.opt_vio(st_idx:ed_idx));
    OPT.time = mean(table.opt_time(st_idx:ed_idx));
    OPTnr.fval = mean(table.opt_wo_fval(st_idx:ed_idx));
    OPTnr.vio = mean(table.opt_wo_vio(st_idx:ed_idx));
    OPTnr.time = mean(table.opt_wo_time(st_idx:ed_idx));
    RDTSH.fval = mean(table.rdtsh_fval(st_idx:ed_idx));
    RDTSH.vio = mean(table.rdtsh_vio(st_idx:ed_idx));
    RDTSH.time = mean(table.rdtsh_time(st_idx:ed_idx));
    TSH.fval = mean(table.tsh_fval(st_idx:ed_idx));
    TSH.vio = mean(table.tsh_vio(st_idx:ed_idx));
    TSH.time = mean(table.tsh_time(st_idx:ed_idx));
    SRIGH.fval = mean(table.srigh_fval(st_idx:ed_idx));
    SRIGH.vio = mean(table.srigh_vio(st_idx:ed_idx));
    SRIGH.time = mean(table.srigh_time(st_idx:ed_idx));
end