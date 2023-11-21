function [X, Y] = HystDataSimulation(N, x_liml, x_limu, phi, s, n)

% Define a range of X_i values
X = linspace(x_liml, x_limu, N); % Adjust the range and number of points as needed
% Set a constant phi value

% Initialize an array to store the function values
N = size(X,2);
Y = zeros(2, N);

% Calculate the function values for each X_i
for i = 1:length(X)
    X_i = X(i);
    
    % Assuming X_i_minus_1 is a constant for this test
    if i == 1
        X_i_minus_1 = X_i-1;
        X_i_plus_1 = X(i + 1);
    elseif i == length(X)
        X_i_minus_1 = X(i - 1);
        X_i_plus_1 = X_i+1;
    else
        X_i_minus_1 = X(i - 1);
        X_i_plus_1 = X(i + 1);
    end
    
    % Calculate the function value using the g function
    Y(1,i) = hysteresis_model(X_i, X_i_minus_1, phi);
    Y(2,i) = hysteresis_model(X_i, X_i_plus_1, phi);
end

% Create a plot

Y(1,:) = smooth(Y(1,:), n,'moving');
Y(2,:) = smooth(Y(2,:), n,'moving');

% Add noise

Y = awgn(Y, s);
end