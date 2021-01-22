%% Calculate battery stage-of-charge.
% Args:
%   y0: total initial capacity in mAh
%   Tc: ambient temperature in Celsius
%   I: average current draw in mA
%   t: target time for SoC evaluation in days
% 
% Return:
%   SoC: estimated SoC

function [SoC] = soc(y0, Tc, I, t)
% required constants
c = 0.56418;
Ea = 1.1949;
A = 0.96397;
R = 0.008314;

% correct initial capacity according to ambient temperature
% valid range for temperature is -5 to 40 Celsius
CF = 1.0;
if Tc >= -5.0 && Tc < 10.0
    CF = (-5.117e-7)*(Tc+5)^3 + 1.0076e-3*(Tc+5) + 9.98e-1;
elseif Tc >= 10.0 && Tc < 25.0
    CF = 2.2375e-6*(Tc-10)^3 + (-2.3027e-5)*(Tc-10)^2 + 6.6620e-4*(Tc-10)^1 + 1.0114;
elseif Tc >= 25.0 && Tc < 32.5
    CF = (-2.0925e-5)*(Tc-25)^3 + 7.7663e-5*(Tc-25)^2 + 1.4817e-3*(Tc-25)^1 + 1.0237;
elseif Tc >= 32.5 && Tc < 40
    CF = 1.7473e-5*(Tc-32.5)^3 + (-3.9315e-4)*(Tc-32.5)^2 + (-8.8444e-4)*(Tc-32.5)^1 + 1.0303;
end
%fprintf("%f\n", CF);
y0 = y0 * CF;
i0 = y0 * c;
j0 = y0 * (1 - c);
k = A * exp(-Ea / (R * (Tc + 273.15))); % note the temperature here is in Kelvin!
%fprintf("%f, %f\n", i0, k);

dt = 1; % in days
tt = t;
% iteratively compute battery lifetime
t0 = 0.0;
while i0 > 0 && t0 < tt
    t0 = t0 + dt;
    y0 = i0 + j0;
    i0 = i0 * exp(-k*dt) + (y0*k*c-I)*(1-exp(-k*dt))/k - I*c*(k*dt-1+exp(-k*dt))/k;
    j0 = j0 * exp(-k*dt) + y0*(1-c)*(1-exp(-k*dt)) - I*(1-c)*(k*dt-1+exp(-k*dt))/k;
    %fprintf("%f, %f\n", t0, i0);
end
%fprintf("%f\n", i0);
SoC = i0 / (y0*c);
end
