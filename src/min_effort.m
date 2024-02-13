function [x, optimal_value] = min_effort(A, y, U, find_x)
    % min_effort Our algorithm for solving the minimum effort problem
    %    minimize     ||x||_infty
    %    subject to   A * x = y1
    %
    % Inputs:
    % A (matrix): matrix A from above.
    % y (vector): vector y from above.
    % U (matrix): set of extremal points of the dual problem computed by get_u(A, B).
    % find_x (function, optional): upon calling find_x(D, d, I0) solves
    %       the minimum effort problem D * x_out = d. Even though
    %       this procedure should handle any case, it may used to handle
    %       problematic cases in a fast way.
    %
    % Outputs:
    % x (vector): optimal solution of the above problem.
    % optimal_value (scalar): optimal value of the above problem.

    % Specify find_x if not provided
    if nargin < 4
        find_x = @(varargin) [];
    end

    tol = 1e-10;
    n = size(A, 2);

    % Compute the optimal dual solution
    [optimal_value, i_max] = max(U*y);
    u_opt = U(i_max, :)';

    % Assign the index sets for complementarity
    I0 = abs(A'*u_opt) <= tol;
    I1 = A'*u_opt > tol;
    I2 = A'*u_opt < -tol;

    % Use the complementarity conditions to compute the primal solution
    x = zeros(n,1);
    x(I1) = optimal_value;
    x(I2) = -optimal_value;

    D = A(:,I0);
    d = y - A(:,I1|I2)*x(I1|I2);

    % We need to solve the following equation for x(I0)
    % D * x(I0) = d;
    x_user = find_x(D, d, I0);
    if ~isequal(x_user, [])
        % Use user-provided solution in present
        x(I0) = x_user;
    elseif cond(D'*D) <= 1e10
        % The simple case when it is well-conditioned
        x(I0) = D \ d;
    else
        % For the general case, we compute the l2 minimal solution.
        % This does not need to be optimal.
        [~, basiccol] = rref(D');
        D = D(basiccol,:);
        d = d(basiccol);
        x(I0) = D' * ((D * D') \ d);
    end

    % Check for solution optimality
    if norm(A*x-y) > tol
        throw("Problem was not solved");
    end
    if max(abs(x)) > optimal_value + tol
        warning("COMPUTATIONS FAILED. SOLUTION IS SOBOPTIMAL.");
    end
end

