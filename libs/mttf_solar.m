%% Computing the MTTF of solar cell
% Args:
%   Ni: the struct at the target grid location
%
% Return:
%   MTTFr: estimated MTTF ratio to standard environment of 25 Celsius

function [MTTFr] = mttf_solar(Ni)
% required constants
Ea = 0.9;           % activation energy in eV
k = 8.617333e-5;    % Boltzmann constant in eV/K
Tref = 25 + 273.15; % reference ambient temperature in K
N_bin = length(Ni.Tcen);  % number of temperature bins

% compute the expectation of MTTFr
MTTFr = 0.0;
for j=1:N_bin
    Tj = Ni.Tcen(j) + 273.15; % convert to Kelvin
    MTTF_cur = exp(Ea / k * (1/Tj - 1/Tref));
    MTTFr = MTTFr + Ni.Tcnt(j) * MTTF_cur;
end
end

