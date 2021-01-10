%% Computing the MTTF of solar cell
% Args:
%   Ni: the struct at the target grid location
%
% Return:
%   MTTFr: estimated MTTF ratio to standard environment of 25 Celsius

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

