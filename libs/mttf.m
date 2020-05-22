function MTTFr = mttf(Tc)
% Calculate mean-time-to-failure of TDDB.
% Args:
%   Tc: core temperature in Kelvin
%
% Return:
%   MTTF: estimated MTTF ratio to 25 Celsius

% required constants
%Ao = 4.0;                            % empirical value
c = 10 * 1.1949e3 / (6.022 * 1.38);  % Ea/k

% referece temperature (25 Celsius) in Kelvin
Tref = 25 + 273.15;

MTTFr = exp(c / Tc) / exp(c / Tref);
%fprintf("%f\n", MTTFr);
end