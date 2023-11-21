function [Ax, bx, Ay, by] = build_constraints(N, lim)
    Ax = [eye(N-1) zeros(N-1,1)] - [zeros(N-1,1) eye(N-1)];
    bx = -lim * ones(N - 1, 1);

    Ay = [Ax zeros(N-1, N)];  % Pad with zeros for the second block
    by = -lim * ones(N - 1, 1);
end