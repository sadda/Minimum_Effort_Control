function [x, optimal_value] = min_effort(A, B, y, U, find_x)
    % min_effort Our algorithm for solving the minimum effort problem
    %    minimize     ||x||_infty
    %    subject to   A * x = y1
    %                 B * x <= y2
    %
    % Inputs:
    % A (matrix): matrix A from above.
    % B (matrix): matrix B from above.
    % y (vector): vector y from above.
    % U (matrix): set of extremal points of the dual problem computed by get_u(A, B).
    % find_x (function, optional): upon calling find_x(A_big, y, I0, J, x) solves
    %       the minimum effort problem A_big(J, I0) * x_out = y(J)
    %       which satisfied all hte conditions.
    %
    % Outputs:
    % x (vector): optimal solution of the above problem.
    % val (scalar): optimal value of the above problem.

    % Specify find_x if not provided
    if nargin < 5
        find_x = @(varargin) [];
    end

    tol = 1e-10;
    A_big = [A; B];
    m1 = size(A, 1);
    m2 = size(B, 1);
    n = size(A_big, 2);

    % Compute the optimal dual solution
    [optimal_value, i_max] = max(U*y);
    u_opt = U(i_max, :)';

    % Assign the index sets for complementarity
    I0 = abs(A_big'*u_opt) <= tol;
    I1 = A_big'*u_opt > tol;
    I2 = A_big'*u_opt < -tol;
    J = [true(m1, 1); u_opt(m1+1:m1+m2) < -tol];

    % Use the complementarity conditions to compute the primal solution
    x = zeros(n,1);
    x(I1) = optimal_value;
    x(I2) = -optimal_value;

    D = A_big(J,I0);
    d = y(J) - A_big(J,I1|I2)*x(I1|I2);

    % We need to solve the following equation for x(I0)
    % D * x(I0) = d;
    x_user = find_x(A_big, y, I0, J, x);
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
    if norm(A*x-y(1:m1)) > tol || (m2 > 0 && norm(max(B*x - y(m1+1:m1+m2), 0)) > tol)
        throw("Problem was not solved");
    end
    if max(abs(x)) > optimal_value + tol
        warning("COMPUTATIONS FAILED. SOLUTION IS SOBOPTIMAL.");
    end
end

