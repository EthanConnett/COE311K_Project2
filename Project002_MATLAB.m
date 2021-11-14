%% COE 311K Project 2


%% Constants
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

%% Problem 1

% % dg / dt = lambdaP * g * (1 - g) - lambdaA * g - lambdaK * g * f
% df / dt = -1 * lambdaD * f + p

% Create the tumor model using the preliminary conditions
delta_vec = [0.001, 0.001, 0.001, 0.001];
tau_vec = [10,20,30,40];
kappa = 3;
params.dt = 0.01*params.T/(2^kappa);
[f, g, p, t_vec] = tum_model(delta_vec,tau_vec, kappa, params);


% Find Errors for given k's, exp in report how low k is no good
errors = zeros(10, 1);
final_g_vals = zeros(10,1);

for kappa = 1:10
    params.dt = 0.01 * params.T / 2^kappa;
    [f, g, p, t_vec] = tum_model(delta_vec,tau_vec, kappa, params);
    
    params.dt = 0.01 * params.T / 2^(kappa-1);
    [f, g2, p, t_vec] = tum_model(delta_vec,tau_vec, kappa-1, params);
    errors(kappa) = 100*((g(end)-g2(end))/g(end));
    final_g_vals(kappa) = g(end);
end
    
%plot(t_vec, f)
%hold
%plot(t_vec, p)
%plot(t_vec, g)

%% Problem 2

% Define Spacing
kappa = 3;

% Define function that needs to be optimized
tum_fxn_sigma = @(x) tum_model(x,tau_vec, kappa, params) ;
J_orig = @(x) j_gen(x, params, tum_fxn_sigma);

% Initial guess for 4 treatments with equal vol fractions
x0 = [0.001, 0.001, 0.001, 0.001];

%Boundary Conditions
x_min = [0,0,0,0];
x_max = [0.01, 0.01, 0.01, 0.01];


% Tolerances
tol_x = 1e-9;
tol_fun = 1e-9;
max_iter = 400;

% Optimize using fmincon
tic
[x_opt, Jval, ~, ~, ~, ~, ~] = fmincon(J_orig, x0, [], [], [], [], ...
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
disp(x_opt)



% Compute solution using intial parameters, plot intial solution
[f_in, g_in, p_in, t_vec_in] = tum_model(delta_vec,tau_vec, kappa, params);

% Compute solution using optimal parameters, plot optimal solution
[f_opt, g_opt, p_opt, t_vec_opt] = tum_model(x_opt,tau_vec, kappa, params);

%Plot results
figure('Name', 'F Function Just Optimizing Vol Frac');
plot(t_vec, f_in)
hold on;
plot(t_vec, f_opt)
legend('Initial Guess', 'Optimal Params');
hold off;

figure('Name', 'G Function Just Optimizing Vol Frac');
plot(t_vec, g_in)
hold on;
plot(t_vec, g_opt)
legend('Initial Guess', 'Optimal Params');
hold off;


%% Problem 3

% Define spacing
kappa = 3;

% Define cost function
% First 4 values are the volume fractions (deltas), next 4 values are the 
% values for time inputs (taus).
tum_fxn_sigma = @(x) tum_model(x(1:4), x(5:8), kappa, params) ;
J_orig = @(x) j_new(x, params, tum_fxn_sigma);

% Initial guess for the four treatments
x0 = [0.001, 0.001, 0.001, 0.001, 9.7955,19.5171,30.6384,40.6257];

% Generate tau constraints
for i = 1:length(tau_vec)
   tau_mins(i) = 5  + 10*(i-1);
   tau_maxes(i) = 15 + 10*(i-1);
end

%Define constraints
x_min = [0,0,0,0, tau_mins];
x_max = [0.01, 0.01, 0.01, 0.01, tau_maxes];

% Define tolerances
tol_x = 1e-9;
tol_fun = 1e-9;
max_iter = 400;

% Optimize using fmincon
tic
[x_opt, Jval, ~, ~, ~, ~, ~] = fmincon(J_orig, x0, [], [], [], [], ...
                                x_min, x_max,[], ...
                                optimset('TolX',tol_x, ...
                                'TolFun', tol_fun, ...
                                'MaxIter', max_iter, ...
                                'Display','iter-detailed'));

toc

% Compute solution using intial parameters, plot intial solution
[f_in, g_in, p_in, t_vec_in] = tum_model(delta_vec,tau_vec, kappa, params);

% Compute solution using optimal parameters, plot optimal solution
[f_opt, g_opt, p_opt, t_vec_opt] = tum_model(x_opt(1:4),x_opt(5:8), kappa, params);

%Plot results
figure('Name', 'F Function Just Optimizing Vol Frac');
plot(t_vec, f_in)
hold on;
plot(t_vec, f_opt)
legend('Initial Guess', 'Optimal Params');
hold off;

figure('Name', 'G Function Just Optimizing Vol Frac');
plot(t_vec, g_in)
hold on;
plot(t_vec, g_opt)
legend('Initial Guess', 'Optimal Params');
hold off;
