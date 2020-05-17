%% test script for mttf and soh bound conversion functions
clc;
clear;
close all;
warning('off','all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test Psoh_bound conversion function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SoHref = 0.8;
T = 3;
Ti = linspace(10, 40, 31);
P_sohi = Psoh_bound(SoHref, T, Ti);
figure;
plot(Ti, P_sohi, 'b-*');
title('Power bound (W) from SoH at various temperature');

Tamb = 35;
T = 2;
nref = 31;
SoHref = linspace(0.6, 0.9, nref);
P_sohr = zeros(nref, 1);
for i=1:nref
    P_sohr(i) = Psoh_bound(SoHref(i), T, Tamb);
end
figure;
plot(SoHref, P_sohr, 'b-*');
title('Power bound (W) from SoH at various SoHref');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test Pmttf_bound conversion function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MTTFref = 0.8;
P_mttfi = Pmttf_bound(MTTFref, Ti);
figure;
plot(Ti, P_mttfi, 'b-*');
title('Power bound (W) from MTTF at various temperature');

Tamb = 30;
nref = 31;
MTTFref = linspace(0.6, 0.9, nref);
P_mttfr = zeros(nref, 1);
for i=1:nref
    P_mttfr(i) = Pmttf_bound(MTTFref(i), Tamb);
end
figure;
plot(MTTFref, P_mttfr, 'b-*');
title('Power bound (W) from MTTF at various MTTFref');