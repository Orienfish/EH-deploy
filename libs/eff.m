%% Compute solar panel's end-to-end conversion efficiency at one location
% Args:
%   Tamb: average ambient temperature
%   A: surface area (m^2) of solar panel
%   dni: average dni at the location
%
% Return:
%   xi: estimated solar panel's conversion efficiency at the location

function [xi] = eff(Tamb, A, dni)
% required constants
a = -6.66e-4;
b = 0.113;
Rth = 0.08;

xi = (a * Tamb + b) / (1 - a * Rth * A * dni);
end

