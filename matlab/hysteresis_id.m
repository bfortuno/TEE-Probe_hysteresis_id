% ------------------------------------------------------------------------ 
% Hysteresis Identification Optimization Script
% ------------------------------------------------------------------------
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

%% Read xlsx
close all
load('Calibration/normal_ap.mat');
load('Calibration/normal_ml.mat');

% Change motor {1/2} to change motor data
T = readtable("Calibration/motor1_current60.xlsx");

current = T.Current;
steps = T.MotorPosition;

tip_q1 = T.TipOrientationX;
tip_q2 = T.TipOrientationY;
tip_q3 = T.TipOrientationZ;
tip_q4 = T.TipOrientationW;

base_q1 = T.BaseOrientationX;
base_q2 = T.BaseOrientationY;
base_q3 = T.BaseOrientationZ;
base_q4 = T.BaseOrientationW;

tip_quat = [tip_q1,tip_q2,tip_q3,tip_q4];
tip_pos = [T.TipPositionX, T.TipPositionY, T.TipPositionZ];
base_quat = [base_q1,base_q2,base_q3,base_q4];
base_pos = [T.BasePositionX, T.BasePositionY, T.BasePositionZ];

tip_PlaneProj = tip_quat(:,1:3)' - normal_ap * (normal_ap' * tip_quat(:,1:3)');
tip_PlaneProj = tip_PlaneProj ./ vecnorm(tip_PlaneProj);

base_PlaneProj = base_quat(:,1:3)' - normal_ap * (normal_ap' * base_quat(:,1:3)');
base_PlaneProj = base_PlaneProj ./ vecnorm(base_PlaneProj);

crossVec = cross(tip_PlaneProj, base_PlaneProj);

theta = sign(normal_ap.' * crossVec) .* atan2d(vecnorm(crossVec), dot(base_PlaneProj, tip_PlaneProj));

%% Current hysteresis plot (only for analysis)
current_smth = 2.69*smooth(current,10);         % 2.69 is the factor of the current
figure(2)
clf
hold on
plot(steps, current_smth)
ylabel("current [mA]")
xlabel("motor position [steps]")
title("ML current hysteresis (moving average window = 10)")
hold off

%% Reorder the date depending on the movement
[Mmax,Imax] = max(steps);
[Mmin,Imin] = min(steps);

n = numel(steps);

% Create steps_pos and steps_neg
allIndices = 1:n;
indicesInRange = mod(Imax:Imin - 1, n) + 1;
indicesOutOfRange = setdiff(allIndices, indicesInRange);

[steps_pos, indicesOutOfRange_2] = sort(steps(indicesOutOfRange));
steps_neg = steps(indicesInRange);

t_pos = theta(indicesOutOfRange);
theta_pos = t_pos(indicesOutOfRange_2);
theta_neg = theta(indicesInRange);

figure(1)
hold on
plot(steps_pos,theta_pos, "Color","r")          % postero-anterior
plot(steps_neg, theta_neg, "Color","b")         % antero-posterior
hold off


%% Read data from txt file
curve = "postero_anterior"; 
% curve = "antero_posterior"; 

file_path = curve + '.txt';

% Use the load function to read data from the text file
data = load(file_path);
data = sortrows(data,1)';

% Assuming your file has two columns, assign them to x and y
Yh = theta_pos;  % Assuming the first column contains x values
Xh = steps_pos';  % Assuming the second column contains y values
X_val = Xh;
Y_val = Yh;
data = [Xh;Yh];
%% X-axis Optimization
n_var = 6;
var_difference = 10;
N = 2*size(Xh,2);

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
xlabel('\theta');
ylabel('Motor Steps');
title('Plot of Motor steps with respect to \theta');
hold off

%% Y-axis Optimization
clc
dx = 10;
dy = 0.1;
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

X_val_ap = data_ap(1,:);
X_val_pa = data_pa(1,:);
Yap = ppval(ss_ap, X_val_ap);
Ypa = ppval(ss_pa, X_val_pa);
Yap_meas = data_ap(2,:);
Ypa_meas = data_pa(2,:);

% Plotting
figure(1)
clf
hold on
grid on;

% Measured data
plot(X_val_ap, Yap_meas, 'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.5); % Orange
plot(X_val_pa, Ypa_meas, 'color', [0, 0.4470, 0.7410], 'LineWidth', 1.5); % Bluish color

% Fitted data
plot(X_val_ap, Yap, 'color', [0.4660, 0.6740, 0.1880], 'LineWidth', 2); % Green
plot(X_val_pa, Ypa, 'color', [0.4940, 0.1840, 0.5560], 'LineWidth', 2); % Purple

% Customize the plot
ylabel('\theta');
xlabel('Motor Steps');
title('Plot of Motor steps with respect to \theta');

% Add legends
legend('Yap_{measured}', 'Ypa_{measured}', 'Yap_{fitted}', 'Ypa_{fitted}', 'Location', 'Best');

% Add annotations or additional information
text(0.5, 1.5, 'Additional Information', 'FontSize', 12, 'Color', 'k');
hold off
