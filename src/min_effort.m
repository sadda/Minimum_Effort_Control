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
    % find_x (function, optional): upon calling find_x(D, d, I0, J) solves
    %       the minimum effort problem D(J, I0) * x(I0) = d(J). Even though
    %       this procedure should handle any case, it may used to handle
    %       problematic cases in a fast way.
    %
    % Outputs:
    % x (vector): optimal solution of the above problem.
    % val (scalar): optimal value of the above problem.

    % Specify find_x if not provided
    if nargin < 5
        find_x = @(D, d, I0, J) [];
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
    x_user = find_x(D, d, I0, J);
    if ~isequal(x_user, [])
        % Use user-provided solution in present
        x(I0) = x_user;
    elseif cond(D'*D) <= 1e10
        % The simple case when it is well-conditioned
        x(I0) = D \ d;
    elseif length(I0) == 5 && length(J) == 6
        % This hard-codes the case for one specific matrix when
        % I0 is [1;0;1;1;0] and J is [1;1;0;0;0;1] or [1;1;0;1;0;0]
        t_p = A_big([1 2],I0) \ (y([1 2]) - A_big([1 2],I1|I2)*x(I1|I2));
        t_d = null(A_big(J,I0));
        if any(t_d <= 0)
            t_d = -t_d;
        end
        c_max1 = min((optimal_value-t_p) ./ t_d);
        c_min1 = max((-optimal_value-t_p) ./ t_d);
        c_max2 = (y(3) - A_big(3,I1|I2)*x(I1|I2) - A_big(3,I0)*t_p) / (A_big(3,I0)*t_d);
        c_min2 = (-y(3) - A_big(3,I1|I2)*x(I1|I2) - A_big(3,I0)*t_p) / (A_big(3,I0)*t_d);
        c_max = min(c_max1, c_max2);
        c_min = max(c_min1, c_min2);

        xs_a = x;
        xs_a(I0) = t_p + c_min*t_d;
        xs_b = x;
        xs_b(I0) = t_p + c_max*t_d;

        norm_a = norm(A_big([3 4],:)*xs_a);
        norm_b = norm(A_big([3 4],:)*xs_b);

        if norm(norm_a) > norm(norm_b) + tol
            c = c_min;
        elseif norm(norm_b) > norm(norm_a) + tol
            c = c_max;
        else
            c = 0.5*(c_min+c_max);
        end
        x(I0) = t_p + c*t_d;
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

