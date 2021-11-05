%% COE 311K Project 2
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
    p = @(t) (i * 1/(sigma * (2 * pi)^(1/2)) * exp(-1 * (abs(tal - t)^2)/(2 * sigma^2)));
end

function solve = discretization(n)

    global p lambdaA lambdaK lambdaD;
     
    h = 1/n;
    
    g = diag(zeros(n));
    h = diag(zeros(n));
    
    for i = 1:n
        g(i+1) = g(i) + h * (lambdaD * g(i) * (1 - g(i)) - lambdaA * g(i) - lambdaK * g(i) * f(i));
        f(i+1) = f(i) + h * (-1 * lambdaD * f(i) + p(i * h));
    end
end