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
n_mttftemp = 41;
Tc = linspace(0, 40, n_mttftemp); % core temperature in Celsius
MTTF = zeros(1, n_mttftemp);
for i = 1:n_mttftemp
    MTTF(i) = mttf(Tc(i) + 273.15); % convert to Kelvin
end
figure;
plot(Tc, MTTF);
title('MTTF ratio under various core temperature');