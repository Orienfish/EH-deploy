%% test script for temperature conversion, MTTF and SoH models 
clc;
clear;
close all;
warning('off','all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test temperature conversion function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ti = linspace(0, 40, 41);
Pi = [0.05, 0.1, 0.15];
Tcore = zeros(length(Ti), length(Pi));
for i=1:length(Ti)
    for j=1:length(Pi)
        Tcore(i, j) = amb2core(Ti(i), Pi(j)) - 273.15;
    end
end
figure;
for j=1:length(Pi) 
    plot(Ti, Tcore(:, j), 'b-*');
    hold on;
end
title('Core temperature (Celsius) at various temperature');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test transmission power consumption under different distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_dist = 41;
dist = linspace(0, 200, n_dist);
Ptx = zeros(1, n_dist);
for i = 1:n_dist
    Ptx(i) = getPtx(dist(i));
end
figure;
plot(dist, Ptx);
title('Average transmission power under various distance');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test MTTF under various core temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_temp = 46;
Tc = linspace(0, 45, n_temp); % core temperature in Celsius
MTTF = zeros(1, n_temp);
MTTF(:) = mttf(Tc + 273.15); % convert to Kelvin
figure('Position', [0 0 250 250]);
plot(Tc, MTTF, '-*', 'LineWidth', 1);
%title('MTTF ratio under various core temperature');
xlabel('Core Temperature (Celsius)'); ylabel('MTTF Ratio');
ax = gca; ax.FontSize = 16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test SoH under various core temperature and target time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use the same core temperature setting as the last test
n_temp = 46;
Tc = linspace(0, 45, n_temp); % core temperature in Celsius
n_time = 1;
T = linspace(5, 5, n_time);
SoH = zeros(n_time, n_temp);
for j = 1:n_time
    SoH(j, :) = soh(Tc + 273.15, T(j)); % convert to Kelvin
end
figure('Position', [0 0 250 250]);
for j = 1:n_time
    plot(Tc, SoH(j, :), 'r-*', 'LineWidth', 1);
    hold on;
end
xlabel('Core Temperature (Celsius)'); ylabel('SoH');
ax = gca; ax.FontSize = 16;
%title('Battery SoH under various core temperature and target time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test Psoh_bound and Pmttf_bound conversion function
% at various temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SoHref = 0.8;
T = 5;
Ti = linspace(0, 40, 41);
P_sohi = Psoh_bound_Tavg(SoHref, T, Ti);
MTTFref = 0.8;
P_mttfi = Pmttf_bound_Tavg(MTTFref, Ti);

figure('Position', [0 0 300 300]);
plot(Ti, P_sohi, '-*', 'LineWidth', 1);
hold on;
plot(Ti, P_mttfi, '-^', 'LineWidth', 1);
xlabel('Ambient Temperature (°C)'); ylabel('Power Bound (W)');
legend({'SoH Constraint', 'MTTF Constraint'}, 'FontSize', 16);
ax = gca; ax.FontSize = 16;
title('Power bound (W) from SoH and MTTF at various temperature');


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
    P_sohr(i) = Psoh_bound_Tavg(SoHref(i), T, Tamb);
end

MTTFref = linspace(0.6, 0.9, nref);
P_mttfr = zeros(nref, 1);
for i=1:nref
    P_mttfr(i) = Pmttf_bound_Tavg(MTTFref(i), Tamb);
end

figure('Position', [0 0 300 300]);
plot(SoHref, P_sohr, '-*', 'LineWidth', 1);
hold on;
plot(MTTFref, P_mttfr, '-^', 'LineWidth', 1);
xlabel('Electronics MTTF or SoH Bounds'); ylabel('Power Bound (W)');
legend({'Battery SoH', 'Electronics MTTF'}, 'FontSize', 16);
ax = gca; ax.FontSize = 16;
grid on;
%title('Power bound (W) from SoH and MTTF at various ref. value');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare single-use batteries and rechargable batteries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nDay = 301;
Day = linspace(0, 600, nDay);
Tamb = 32; pwr = 0.1;
Tcell = amb2core(Tamb, pwr);
y0 = 1000;  % 1000 mAh
I = 10;     % mA
bsoh = zeros(nDay, 1);
bsoc = zeros(nDay, 1);
for i=1:nDay
    bsoc(i) = soc(y0, Tamb, I, Day(i));
    bsoh(i) = soh(Tcell, Day(i) / 365);
end
figure('Position', [0 0 400 300]);
plot(Day, bsoc, '-s', 'Linewidth', 1);
hold on;
plot(Day, bsoh, '-*', 'Linewidth', 1);
legend({'SoC of single-use battery', 'SoH of rechargeable battery'}, 'FontSize', 16, 'Location', 'southeast');
ylim([0, 1.1]); xlabel('Elapsed Time (day)'); ylabel('Battery Status');
ax = gca; ax.FontSize = 16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test solar panel conversion efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_temp = 51;
Tamb = linspace(0, 50, n_temp); % ambient temperature in Celsius
xi = zeros(n_temp, 1);
for i=1:n_temp
    xi(i) = eff(Tamb(i), 0.01, 300);
end
figure;
plot(Tamb, xi, '-s', 'Linewidth', 1);
xlabel('Average ambient temperature (Celsius)'); ylabel('Conversion efficiency');
ax = gca; ax.FontSize = 16;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test solar panel MTTF in ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_temp = 31;
Tamb = linspace(20, 50, n_temp); % ambient temperature in Celsius
MTTFr = zeros(n_temp, 1);
for i=1:n_temp
    MTTFr(i) = mttf_solar(Tamb(i));
end
figure;
plot(Tamb, MTTFr, '-s', 'Linewidth', 1);
xlabel('Average ambient temperature (Celsius)'); ylabel('Solar panel MTTF in ratio');
ax = gca; ax.FontSize = 16;

% This is the function using single ambient temperature as input
% while 
function [MTTFr] = mttf_solar(Tamb)
% required constants
Ea = 0.9;           % activation energy in eV
k = 8.617333e-5;    % Boltzmann constant in eV/K
Tref = 25 + 273.15; % reference ambient temperature in K
%N_bin = length(Ni.Tcen);  % number of temperature bins

Tamb = Tamb + 273.15;
% Compute the expectation of MTTFr
MTTFr = exp(Ea / k * (1/Tamb - 1/Tref));
end