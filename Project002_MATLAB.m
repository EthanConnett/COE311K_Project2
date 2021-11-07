%% COE 311K Project 2
clear all
% dg / dt = lambdaP * g * (1 - g) - lambdaA * g - lambdaK * g * f
% df / dt = -1 * lambdaD * f + p

%% Constants
global p lambdaP lambdaA lambdaK lambdaD T g0 f0 NTreat sigma a b;

T = 50;
lambdaP = 0.1;
lambdaA = 0.005;
lambdaK = 4;
lambdaD = 0.05;
g0 = 0.01;
f0 = 0.001;
NTreat = 4;
sigma = 2;
a = 0.0005;
b = 0.00005;

for i = 1:NTreat
    tal = i * 10;
    p_old = @(t) 0;
    p = @(t) p_old(t) + (i * 1/(sigma * (2 * pi)^(1/2)) * exp((-1 * (abs(tal - t)^2)/(2 * sigma^2))));
    p_old = p;
end

%ek = errorG();

% for k = 0:10
%     dt = 0.01 * T / 2^k;
%     [f, g] = discretization(dt, g0, f0);
% end