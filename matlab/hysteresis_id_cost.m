function    [cost,Ys, ss]    =   hysteresis_id_cost(x, Xh, Yh, X_val, Y_val, solver, type)
% hysteresis_id_cost - Computation of squared error between model and
% measurements
%
% Syntax:  [cost,Ys, ss] = kinetostatic_id_cost(x, Xh, Yh, X_val, Y_val, solver, type)
%
% Inputs:
%    x              -   vector corresponding to the index in X to create model to compare,
%                       the form is [1x2n]
%    Xh, Yh         -   Values of train set
%    X_val, Y_val   -   Values of Validation set
%    solver         -   Can be
%   
% Outputs:
%    output1 - Description
%    output2 - Description
%
% Example: 
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% Author: Benjamin Fortuno
% NearLab, Politecnico di Milano
% email address: benjaminignacio.fortuno@mail.polimi.it
% November 2023; Last revision: 17-Nov-2023

N = size(x,2);

if strcmp(type,'XF')
    N = 2 * N;
    Yf = Yh(:, x(:,1:N/2));
elseif strcmp(type,'YF')
    Yf = x(:,N/2+1:end);
end

if strcmp(solver,'splines')
    Ys = spline(Xh(:, x(:,1:N/2)), Yf, X_val);
    ss = spline(Xh(:, x(:,1:N/2)), Yf);
elseif strcmp(solver,'makima')
    Ys = makima(Xh(:, x(:,1:N/2)), Yf, X_val);
    ss = makima(Xh(:, x(:,1:N/2)), Yf);
elseif strcmp(solver,'pchip')
    Ys = pchip(Xh(:, x(:,1:N/2)), Yf, X_val);
    ss = pchip(Xh(:, x(:,1:N/2)), Yf);
end

err_vec = Yh - Ys;
% assignin("base", "x", x);

%% Plot
figure(1);
clf
hold on
grid on;
c1 = plot(X_val, Y_val, 'r', 'LineWidth', 1.5);
c2 = plot(X_val,Ys,'b','LineWidth',2);
xlabel('\theta_i');
ylabel('T(\theta_i)');
title('Plot of T(\theta_i) for a Constant \phi');
xlim([min(X_val) - 5, max(X_val) + 5]);
ylim([min(Y_val) - 0.02, max(Y_val) + 0.02]);
hold off

%% Compute sum of squared errors
cost    =   sum(err_vec.*err_vec);

end