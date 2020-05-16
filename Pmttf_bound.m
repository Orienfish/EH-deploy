%% Convert the lower bound on SoH to the upper bound on power 
% Args:
%   MTTFref: lower bound target for MTTF, a ratio
%   Ti: list of average ambient temperatures in Celsius at the grid locations
%
% Return:
%   P_mttfi: list of upper bound for average power

function P_mttfi = Pmttf_bound(MTTFref, Ti)
% coefficients for temperature conversion
% Tcore (Kelvin) = k_1 * pwr (W) + k_2 * Tamb (Kelvin) + k_3;
k_1 = 4.43815458; % values need to be changed!
k_2 = 1.41558165;
k_3 = 1.7971150984226512;

% referece temperature (25 Celsius) in Kelvin
Tref = 25 + 273.15;

% in the MTTF model
c = 10 * 1.1949e3 / (6.022 * 1.38);

P_mttfi = (1/k_1) * ((Tref*c)/(Tref*log(MTTFref)+c) - k_2 Ti - k_3);
end

