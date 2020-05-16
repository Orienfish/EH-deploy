%% Convert the lower bound on SoH to the upper bound on power 
% Args:
%   SoHref: lower bound target for SoH, 0 <= SoHref <= 1
%   T: the length of the time period in years
%   Ti: list of average ambient temperatures at the grid locations
%
% Return:
%   P_soh: the upper bound for average power

function P_soh = Psoh_bound(SoHref, T, Ti)
% coefficients for temperature conversion
% Tcore (Celsius) = k_1 * pwr (W) + k_2 * Tamb (Celsius) + k_3;
k_1 = 4.43815458;
k_2 = 1.41558165;
k_3 = 1.7971150984226512;

% referece temperature (25 Celsius) in Kelvin
Tref = 25 + 273.15;

% convert target time from years to seconds
Tsec = T * 365 * 24 * 3600;

% coefficients in the battery SoH model
kt = 4.14e-10;          % time stress coefficient in s-1
kT = 6.93e-2;           % temperature stress coefficient

denom = 1 - (1/(kT*Tref)) * log(-log(SoHref)/(kt*Tsec));
P_soh = (1/k_1) * ((Tref / denom) - k_2 * Ti - k_3);
end

