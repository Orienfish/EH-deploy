%% Convert ambient temperature to core temperature
% Args:
%   Ti: average ambient temperatures in Celsius at the grid location
%   Pi: average power in W at the grid location
%
% Return:
%   Tcorei: average core temperature in Kelvin at the grid location

function Tcorei = amb2core(Ti, Pi)
% coefficients for temperature conversion
global k_1;
global k_2;
global k_3;
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

% convert ambient temperature from Celsius to Kelvin
Ti = Ti + 273.15;

% calculate the core temperature in Kelvin
Tcorei = k_1 * Pi + k_2 * Ti + k_3;
end

