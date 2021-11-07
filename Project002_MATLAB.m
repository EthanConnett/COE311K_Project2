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

%% 1 Setup

for i = 1:NTreat
    tal = i * 10;
    p_old = @(t) 0;
    p = @(t) p_old(t) + (i * 1/(sigma * (2 * pi)^(1/2)) * exp((-1 * (abs(tal - t)^2)/(2 * sigma^2))));
    p_old = p;
end

<<<<<<< HEAD
for k = 0:10
    dt = 0.01 * T / 2^k;
    [f, g] = discretization(dt, g0, f0);
    
end

%% Problem 1

% Discretize

% Find Errors for given k's, exp in report how low k is no good

% Fix k such that error is less than 1

%ek = errorG();

% for k = 0:10
%     dt = 0.01 * T / 2^k;
%     [f, g] = discretization(dt, g0, f0);
% end

%% Problem 2

% what is x again:
t = 0;
x = 0;

g = @(t) t;
J_orig = @(x) j_gen(x, a, b, g);

% Initial guess
x0 = [0.001, 0.001, 0.001, 0.001];

x_min = x(1);
x_max = x(end);

% Prof said not to change this
tol_x = 1e-9;
tol_fun = 1e-9;
max_iter = 400;

tic
[xopt, Jval, ~, ~, ~, ~, ~] = fmincon(J_orig, x0, [], [], [], [], ...
                                x_min, x_max,[], ...
                                optimset('TolX',tol_x, ...
                                'TolFun', tol_fun, ...
                                'MaxIter', max_iter, ...
                                'Display','iter-detailed'));

toc

% Display results
disp('initial guess')
disp(x0)
disp('optimal parameters')

% Compute solution using intial parameters
[g0, f0, tum_init] = F(x0);

% Compute solution using intial parameters
[gopt, fopt, tum_opt] = F(x0);

% Plot
plot(tum_int);
hold on;
plot(tim_out);
legend("Intial guess", "Optimal");

