%% test script for mttf and soh bound conversion functions
clc;
clear;
close all;
warning('off','all');
addpath('./libs');

% call the amb2core function to load the global variables
Tcorei = amb2core(25, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test Psoh_bound and Pmttf_bound conversion function
% at various temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SoHref = 0.8;
T = 5;
Ti = linspace(0, 40, 41);
P_sohi = Psoh_bound(SoHref, T, Ti);
figure('Position', [0 0 300 300]);
plot(Ti, P_sohi, '-*');
xlabel('Ambient Temperature (°C)'); ylabel('Power Bound (W)');
%title('Power bound (W) various temperature');
hold on;
MTTFref = 0.8;
P_mttfi = Pmttf_bound(MTTFref, Ti);
%figure;
plot(Ti, P_mttfi, '-^');
legend({'SoH Constraint', 'MTTF Constraint'}, 'FontSize', 16);
ax = gca; ax.FontSize = 16;
%title('Power bound (W) from MTTF at various temperature');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test Psoh_bound and Pmttf_bound conversion function
% at various reference values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tamb = 35;
T = 5;
nref = 31;
SoHref = linspace(0.6, 0.9, nref);
P_sohr = zeros(nref, 1);
for i=1:nref
    P_sohr(i) = Psoh_bound(SoHref(i), T, Tamb);
end
figure('Position', [0 0 300 300]);
plot(SoHref, P_sohr, '-*');
%title('Power bound (W) from SoH at various SoHref');
xlabel('Target'); ylabel('Power Bound (W)');
hold on;
MTTFref = linspace(0.6, 0.9, nref);
P_mttfr = zeros(nref, 1);
for i=1:nref
    P_mttfr(i) = Pmttf_bound(MTTFref(i), Tamb);
end
%figure;
plot(MTTFref, P_mttfr, '-^');
legend({'SoH Constraint', 'MTTF Constraint'}, 'FontSize', 16);
ax = gca; ax.FontSize = 16;
%title('Power bound (W) from MTTF at various MTTFref');