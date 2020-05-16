%% test script for mttf and soh bound conversion functions
clc;
clear;
close all;
warning('off','all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test Psoh_bound conversion function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SoHref = 0.8;
T = 3;
Ti = linspace(10, 40, 31);
P_sohi = Psoh_bound(SoHref, T, Ti);
figure;
plot(Ti, P_sohi);
title('Power bound (W) from SoH constraint');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test Pmttf_bound conversion function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MTTFref = 0.8;
P_mttfi = Pmttf_bound(MTTFref, Ti);
figure;
plot(Ti, P_mttfi);
title('Power bound (W) from MTTF constraint')