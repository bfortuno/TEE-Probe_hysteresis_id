function [var, exitflag] = BP_solver(data_ap, data_ml, var0)
% BP_SOLVER solves an optimization problem to estimate parameters.
%
% [VAR, EXITFLAG] = BP_SOLVER(DATA_AP, DATA_ML, VAR0) takes three input
% arguments:
%   - DATA_AP: Data for the first set of planes.
%   - DATA_ML: Data for the second set of planes.
%   - VAR0: Initial guess for the optimization.
%
% The function returns:
%   - VAR: Optimized parameters.
%   - EXITFLAG: Exit flag from the optimization algorithm.

    % Set optimization options
    myoptions = myoptimset;
    myoptions.Hessmethod = 'BFGS';
    myoptions.gradmethod = 'CD';
    myoptions.graddx = 2^-17;
    myoptions.tolgrad = 1e-20;
    myoptions.tolx = 1e-40;
    myoptions.tolfun = 1e-40;
    myoptions.ls_tkmax = 1;
    myoptions.ls_beta = 0.9;
    myoptions.ls_c = 1e-1;
    myoptions.ls_nitermax = 1000;

    % Run solver
    tic
    [var, ~, ~, exitflag] = myfminunc(@(var)optim_f(var), var0, myoptions);
    toc

    %% Solver Function
    function [res] = optim_f(var)
        % Extract normals and error vectors for the two sets of planes
        normal_ap = var(1:3, 1);
        normal_ap = normal_ap / norm(normal_ap);
        err_vec_ap = transpose(normal_ap) * data_ap;

        normal_ml = var(4:6, 1);
        normal_ml = normal_ml / norm(normal_ml);
        err_vec_ml = transpose(normal_ml) * data_ml;

        % Define the objective function (sum of squared errors and
        % perpendicular planes)
        res = sum(err_vec_ap .* err_vec_ap) + sum(err_vec_ml .* err_vec_ml);
        res = res + 6500 * dot(normal_ap, normal_ml);

        %% Error Definition for Evaluation
        % Compute mean and maximum errors for each set of planes
        err_ap = mean(abs(err_vec_ap));
        err_ml = mean(abs(err_vec_ml));
        max_ap_err = max(abs(err_vec_ap));
        max_ml_err = max(abs(err_vec_ml));

        % Store errors in the base workspace for evaluation
        assignin('base', 'err_ap', err_ap);
        assignin('base', 'err_ml', err_ml);
        assignin('base', 'max_ap_err', max_ap_err);
        assignin('base', 'max_ml_err', max_ml_err);
    end
end
