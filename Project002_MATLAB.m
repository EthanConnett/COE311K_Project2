%% COE 311K Project 2
clear all
% dg / dt = lambdaP * g * (1 - g) - lambdaA * g - lambdaK * g * f
% df / dt = -1 * lambdaD * f + p

%% Constants
global lambdaP lambdaA lambdaK lambdaD T g0 f0 NTreat sigma a b;

params.T = 50;
params.lambdaP = 0.1;
params.lambdaA = 0.005;
params.lambdaK = 4;
params.lambdaD = 0.05;
params.g0 = 0.01;
params.f0 = 0.001;
params.NTreat = 4;
params.sigma = 2;
params.a = 0.0005;
params.b = 0.00005;
params.c = 0.001;

%% 1 Setup

% p_old = @(t) 0;
% 
% for i = 1:NTreat
%     tal = i * 10;
%     p = @(t) p_old(t) + (i * 1/(sigma * (2 * pi)^(1/2)) * exp((-1 * (abs(tal - t)^2)/(2 * sigma^2))));
%     p_old = p;
% end
% 
% 
% for k = 0:10
%     dt = 0.01 * T / 2^k;
%     [f, g] = discretization(dt, g0, f0, p);
%     
% end
% 
% 
% dt = 0.01 * T / 2^1;
% [f_tes, g_tes] = discretization(dt, g0, f0, p);
    


%% Problem 1

% Discretize
delta_vec = [0.001, 0.001, 0.001, 0.001];
tau_vec = [10,20,30,40];
kappa = 3;
params.dt = 0.01*params.T/(2^kappa);
[f, g, p, t_vec] = tum_model(delta_vec,tau_vec, kappa, params);


% Find Errors for given k's, exp in report how low k is no good

% Fix k such that error is less than 1

errors = zeros(10, 1);

for kappa = 1:10
    params.dt = 0.01 * params.T / 2^kappa;
    [f, g, p, t_vec] = tum_model(delta_vec,tau_vec, kappa, params);
    
    params.dt = 0.01 * params.T / 2^(kappa-1);
    [f, g2, p, t_vec] = tum_model(delta_vec,tau_vec, kappa-1, params);
    errors(kappa) = 100*((g(end)-g2(end))/g(end));
end
    
%plot(t_vec, f)
%hold
%plot(t_vec, p)
%plot(t_vec, g)

%% Problem 2

% what is x again:

kappa = 3;

tum_fxn_sigma = @(x) tum_model(x,tau_vec, kappa, params) ;
J_orig = @(x) j_gen(x, params, tum_fxn_sigma);

% Initial guess for 4 treatments with equal vol fractions
x0 = [0.001, 0.001, 0.001, 0.001];

x_min = [0,0,0,0];
x_max = [0.01, 0.01, 0.01, 0.01];

j_test = J_orig(x0);

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
disp(xopt)



% Compute solution using intial parameters
% [g0, f0, tum_init] = F(x0);

% Compute solution using intial parameters
% [gopt, fopt, tum_opt] = F(x0);

%

% Plot
% plot(tum_int);
% hold on;
% plot(tim_out);
% legend("Intial guess", "Optimal");

%% Problem 3

% what is x again:

kappa = 3;

tum_fxn_sigma = @(x) tum_model( x(1, :), x(2,:), kappa, params) ;
J_orig = @(x) j_new(x, params, tum_fxn_sigma);

% Initial guess for 4 treatments with equal vol fractions
x0 = [0.001, 0.001, 0.001, 0.001;
      9.7955,19.5171,30.6384,40.6257];

for i = 1:length(tau_vec)
   tau_mins(i) = 5  + 10*(i-1);
   tau_maxes(i) = 15 + 10*(i-1);
end
x_min = [0,0,0,0; tau_mins];
x_max = [0.01, 0.01, 0.01, 0.01; tau_maxes];

j_test = J_orig(x0);

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
disp(xopt)
