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

%% Read data from txt file
curve = "postero_anterior"; 
% curve = "antero_posterior"; 

file_path = curve + '.txt';

% Use the load function to read data from the text file
data = load(file_path);
data = sortrows(data,1)';

% Assuming your file has two columns, assign them to x and y
Xh = data(1, :);  % Assuming the first column contains x values
Yh = data(2, :);  % Assuming the second column contains y values
X_val = Xh;
Y_val = Yh;

%% X-axis Optimization
n_var = 8;
var_difference = 3;

% Choose interpolation method
solver = 'splines'; % You can change this value
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
xlabel('\theta');
ylabel('Motor Steps');
title('Plot of Motor steps with respect to \theta');
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
xlabel('\theta');
ylabel('Motor Steps');
title('Plot of Motor steps with respect to \theta');
hold off

%% Save optimized parameters in 'PPModel.mat'
% Navigate to the parent directory (assuming the "python" folder is at the same level)
parentDir = fileparts(pwd);  % Get the current directory and navigate to its parent
pythonFolderPath = fullfile(parentDir, 'python');

% Save the variable in the "python" folder
file_path = curve + '_PPModel.mat';
save(fullfile(pythonFolderPath, file_path), 'ss');
save(file_path, 'ss');
file_path = curve + '_measurements.mat';
save(file_path, 'data');

%% Plot everything together
ss_ap = load("antero_posterior_PPModel.mat").ss;
ss_pa = load("postero_anterior_PPModel.mat").ss;
data_ap = load("antero_posterior_measurements.mat").data;
data_pa = load("postero_anterior_measurements.mat").data;


Yap = ppval(ss_ap, X_val);
Ypa = ppval(ss_pa, X_val);
Yap_meas = data_ap(2,:);
Ypa_meas = data_pa(2,:);

figure(1)
clf
hold on
grid on;
plot(X_val, Yap_meas, 'r', 'LineWidth', 1.5);
plot(X_val, Yap, 'b', 'LineWidth', 2);
plot(X_val, Ypa_meas, 'r', 'LineWidth', 1.5);
plot(X_val, Ypa, 'b', 'LineWidth', 2);
xlabel('\theta');
ylabel('Motor Steps');
title('Plot of Motor steps with respect to \theta');
hold off
