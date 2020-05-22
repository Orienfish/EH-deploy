%% Calculate battery state-of-health.
% Args:
%   Tcell: battery cell temperature in Kelvin
%   T: the target time in years
%
% Return:
%   SoH: estimated SoH at reference temperature of 25 Celsius

function SoH = soh(Tcell, T)
% referece temperature (25 Celsius) in Kelvin
Tref = 25 + 273.15;

% convert target time from years to seconds
Tsec = T * 365 * 24 * 3600;

% coefficients in the battery SoH model
kt = 4.14e-10;           % time stress coefficient in s-1
kT = 6.93e-2;            % temperature stress coefficient

SoH = exp(-kt * Tsec * exp(kT * Tref * (1 - (Tref / Tcell))));
%fprintf('%f\n', SoH);
end