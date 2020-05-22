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
% test MTTF under various core temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_temp = 51;
Tc = linspace(0, 50, n_temp); % core temperature in Celsius
MTTF = zeros(1, n_temp);
MTTF(:) = mttf(Tc + 273.15); % convert to Kelvin
figure;
plot(Tc, MTTF);
title('MTTF ratio under various core temperature');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test SoH under various core temperature and target time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use the same core temperature setting as the last test
n_time = 10;
T = linspace(1, 5, n_time);
SoH = zeros(n_time, n_temp);
for j = 1:n_time
    SoH(j, :) = soh(Tc + 273.15, T(j)); % convert to Kelvin
end
figure;
for j = 1:n_time
    plot(Tc, SoH);
    hold on;
end
title('Battery SoH under various core temperature and target time');