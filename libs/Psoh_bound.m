%% Convert the lower bound on SoH to the upper bound on power 
% Args:
%   SoHref: lower bound target for SoH at T, 0 <= SoHref <= 1
%   T: the length of the time period in years
%   Ti: list of average ambient temperatures in Celsius at the grid locations
%
% Return:
%   P_sohi: list of upper bound for average power

function P_sohi = Psoh_bound(SoHref, T, Ti)
% coefficients for temperature conversion
% Tcore (Kelvin) = k_1 * pwr (W) + k_2 * Tamb (Kelvin) + k_3;
addpath('./libs');
global k_1;
global k_2;
global k_3;

% referece temperature (25 Celsius) in Kelvin
Tref = 25 + 273.15;

% convert ambient temperature from Celsius to Kelvin
Ti = Ti + 273.15;

% convert target time from years to seconds
Tsec = T * 365 * 24 * 3600;

% coefficients in the battery SoH model
kt = 4.14e-10;           % time stress coefficient in s-1
kT = 6.93e-2;            % temperature stress coefficient

T_celli = Tref / (1 - (1/(kT*Tref)) * log(-log(SoHref)/(kt*Tsec)));
%fprintf('%f %f\n', T_celli, T_celli-273.15);
P_sohi = (1/k_1) * (T_celli - k_2*Ti - k_3);
P_sohi = max(P_sohi, 0); % validity check
end

