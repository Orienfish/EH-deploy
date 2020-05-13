%% Get the power inequality constraint to satisfy ENO
% Args:
%   eta: the sampling frequency
%   B: bandwidth B/s
%   v_cnt: the total number of variables
%   xform: the format of the variables
%   N: the struct of grid locations
%   c: the location of the sink
%
% Reture:
%   A, b: a N_cnt * v_cnt matrix and a N_cnt * 1 vector for power ENO
%         constraint

function [A, b] = power_eno(eta, B, v_cnt, xform, N, c)
% initialization
N_cnt = xform.x_cnt;         % number of grid locations
dist = getDist(N);           % get the distance matrix

% constant parameters setting
P0 = 0.01;                  % sleep power (W)
Es = 0.2*0.2;               % sensing power (W*s)
Prx = 0.1;                  % reception power (W)

% get the transformation matrix to extract s
T_s = [zeros(xform.s_cnt, xform.x_cnt), eye(xform.s_cnt), ...
    zeros(xform.s_cnt, v_cnt-xform.s_end)];
% N_cnt * v_cnt matrix for P0 + si*eta*Es
A = repmat(P0, N_cnt, v_cnt) + eta * Es * T_s;

% set the coefficients for the flow part
for i = 1:N_cnt
    for j = 1:N_cnt
        % outgoing flow from node i to node j, fij
        fij_idx = xform.fij_base + (i-1)*N_cnt + j;
        A(i, fij_idx) = (1/B) * getPtx(dist(i, j));
        % incoming flow from node j to node i, fji
        fji_idx = xform.fij_base + (j-1)*N_cnt + i;
        A(i, fji_idx) = (1/B) * Prx;
    end
    % outgoing flow from node i to the sink, fiB
    fiB_idx = xform.fiB_base + i;
    A(i, fiB_idx) = (1/B) * getPtx(norm(N(i).position - c));
end
fprintf('power eno:\n');
disp(A);
b = vertcat(N(:).Ri);
end

%% Get the transmission power from distance
% Args:
%   dist_m: the distance in m
%
% Return:
%   Ptx: the transmission power in W

function Ptx = getPtx(dist_m)
Pto = 0.12;                 % transmission power baseline
alpha = 3.5;                % distance coefficient
beta = 1e-8;                % distance coefficient

Ptx = Pto + beta * dist_m ^ alpha;
end

%% Get distance matrix between grid points
% Args:
%   N: the struct of grid locations
%
% Return:
%   dist: the N_cnt * N_cnt distance matrix

function dist = getDist(N)
N_cnt = size(N, 1);
dist = zeros(N_cnt);
for i = 1:N_cnt
    for j = i+1:N_cnt
        d = norm(N(i).position - N(j).position);
        dist(i, j) = d;
        dist(j, i) = d;
    end
end
end


