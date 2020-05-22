%% Calculate mean-time-to-failure in a ratio to standard case (25 Celsius).
% Args:
%   Tc: core temperature in Kelvin
%
% Return:
%   MTTF: estimated MTTF ratio to standard environment of 25 Celsius

function MTTFr = mttf(Tc)
% required constants
%Ao = 4.0;                            % empirical value
c = 10 * 1.1949e3 / (6.022 * 1.38);  % Ea/k

% referece temperature (25 Celsius) in Kelvin
Tref = 25 + 273.15;

MTTFr = exp(c / Tc) / exp(c / Tref);
%fprintf("%f\n", MTTFr);
end