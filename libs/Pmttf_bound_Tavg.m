%% Convert the lower bound on SoH to the upper bound on power using average
%  of ambient temperature.
% Args:
%   MTTFref: lower bound target for MTTF, a ratio
%   Ti: list of average ambient temperatures in Celsius at the grid locations
%
% Return:
%   P_mttfi: list of upper bound for average power

function P_mttfi = Pmttf_bound_Tavg(MTTFref, Ti)
% coefficients for temperature conversion
% Tcore (Kelvin) = k_1 * pwr (W) + k_2 * Tamb (Kelvin) + k_3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temperature model for RPi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%k_1 = 4.43815458;
%k_2 = 1.41558165;
%k_3 = -111.71901382193846;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temperature model for low-power devices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_1 = 100;
k_2 = 1.0;
k_3 = 0.0;

% referece temperature (25 Celsius) in Kelvin
Tref = 25 + 273.15;

% convert ambient temperature from Celsius to Kelvin
Ti = Ti + 273.15;

% Ea/k in the MTTF model
c = 10 * 1.1949e3 / (6.022 * 1.38);

T_celli = (Tref * c) / (Tref * log(MTTFref) + c);
%fprintf('%f %f\n', T_celli, T_celli-273.15);
P_mttfi = (1/k_1) * (T_celli - k_2 * Ti - k_3);
P_mttfi = max(P_mttfi, 0);
end
