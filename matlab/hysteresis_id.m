% Hysteresis Identification Optimization Script
%
% Author: Benjamin Fortuno
% Contact: benjaminignacio.fortuno@mail.polimi
%
% Description:
% This script performs hysteresis identification optimization on a
% simulated measurement model. It optimizes the X-axis and Y-axis for a
% given range of theta_i values using surrogate optimization.
%
% Instructions:
% 1. Run the script in MATLAB to perform the optimization.
% 2. Adjust parameters such as phi, solver, n_var, var_difference, dx, dy,
%    and other optimization settings as needed.
%
% Dependencies:
% - Functions: HystDataSimulation, hysteresis_id_cost, surrogateopt
% - MATLAB Optimization Toolbox
%
% Note: Ensure that the necessary functions and toolbox are available in
% the MATLAB environment.

% Clear workspace, command window, and close all figures
clc;
clear all;
close all;

%% Measurements model
N = 400;

% Define a range of theta_i values
x_liml = -55;
x_limu = 99;
phi = pi/3; % You can change this value
s = 55;
n = 10;

% Generate KinetoDataSimulation
[theta_i_values, g_values] = HystDataSimulation(N, x_liml, x_limu, phi, s, n);

% Separate theta and g values for X and Y axis optimization
Xh = theta_i_values(1:2:end);
X_val = theta_i_values(2:2:end);
Yh = g_values(1, 1:2:end);
Y_val = g_values(1, 2:2:end);

%% X-axis Optimization
n_var = 8;
var_difference = 3;

% Choose interpolation method
solver = 'makima'; % You can change this value
x0 = round((0:n_var-1) * (N / 2 - 1) / (n_var - 1)) + 1;
type = 'XF'; % X FIT

% Build constraints
[Ax, bx, Ay, by] = build_constraints(n_var, var_difference);
lb = ones(1, n_var);
ub = ones(1, n_var) * N/2;

% Optimize X-axis
[cost, Ys] = hysteresis_id_cost(x0, Xh, Yh, X_val, Y_val, solver, type);
rng default % For reproducibility
options = optimoptions('surrogateopt');
options.InitialPoints = x0;
options.MaxFunctionEvaluations = 1000; 

objconstr = @(x) hysteresis_id_cost(x, Xh, Yh, X_val, Y_val, solver, type);
x_star_1 = surrogateopt(objconstr, lb, ub, 1:n_var, Ax, bx, [], [], options);
[cost,Ys, ss]    =   hysteresis_id_cost(x_star_1, Xh, Yh, X_val, Y_val, solver, type);

% Plot results
figure(1)
clf
hold on
grid on;
plot(X_val, Y_val, 'r', 'LineWidth', 1.5);
plot(X_val, Ys, 'b', 'LineWidth', 2);
xlabel('\theta_i');
ylabel('T(\theta_i)');
title('Plot of T(\theta_i) for a Constant \phi');
hold off

%% Y-axis Optimization
clc
dx = 7;
dy = 0.02;
yl = Yh(:, x_star_1) - dy;
yu = Yh(:, x_star_1) + dy;
Xl = max(1, x_star_1 - dx);
Xu = min(N/2, x_star_1 + dx);

x0 = [x_star_1 Yh(:, x_star_1)];
type = 'YF'; % Y FIT

% Optimize Y-axis
rng default % For reproducibility
options = optimoptions('surrogateopt');
options.InitialPoints = x0;
options.MaxFunctionEvaluations = 1000; 

objconstr = @(x) hysteresis_id_cost(x, Xh, Yh, X_val, Y_val, solver, type);
lb = [Xl yl];
ub = [Xu yu];

x_star_2 = surrogateopt(objconstr, lb, ub, 1:n_var, Ay, by, [], [], options);
[cost,Ys, ss]    =   hysteresis_id_cost(x_star_2, Xh, Yh, X_val, Y_val, solver, type);

% Plot results
figure(1)
clf
hold on
grid on;
plot(X_val, Y_val, 'r', 'LineWidth', 1.5);
plot(X_val, Ys, 'b', 'LineWidth', 2);
xlabel('\theta_i');
ylabel('T(\theta_i)');
title('Plot of T(\theta_i) for a Constant \phi');
hold off

%% Save optimized parameters in 'PPModel.mat'
% Navigate to the parent directory (assuming the "python" folder is at the same level)
parentDir = fileparts(pwd);  % Get the current directory and navigate to its parent
pythonFolderPath = fullfile(parentDir, 'python');

% Save the variable in the "python" folder
save(fullfile(pythonFolderPath, 'PPModel.mat'), 'ss');
