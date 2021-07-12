function [x, val] = min_effort_linprog(A, y)
    % min_effort_linprog Slow method for solving the minimum effort problem
    %    minimize     ||x||_infty
    %    subject to   Ax = y
    % Use min_effort instead.
    %
    % Inputs:
    % A (matrix): matrix A from above.
    % y (vector): vector y from above.
    % U (matrix): set of extremal points of the dual problem computed by get_u(A).
    %
    % Outputs:
    % x (vector): optimal solution of the above problem.
    % val (scalar): optimal value of the above problem.
    
    [m, n] = size(A);
    
    c      = zeros(n+1,1);
    c(end) = 1;
    A_ineq = [eye(n), -ones(n,1); -eye(n), -ones(n,1)];
    b_ineq = [zeros(n,1); zeros(n,1)];
    A_eq   = [A, zeros(m,1)];
    b_eq   = y;
    lb     = [];
    ub     = [];
    opts   = optimoptions('linprog', 'Display', 'off');
    
    [xz, val] = linprog(c, A_ineq, b_ineq, A_eq, b_eq, lb, ub, opts);
    x = xz(1:n);
end

