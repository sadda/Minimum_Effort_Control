function [x, val] = min_effort(A1, A2, y, U)
    % min_effort Our algorithm for solving the minimum effort problem
    %    minimize     ||x||_infty
    %    subject to   Ax = y
    %
    % Inputs:
    % A (matrix): matrix A from above.
    % y (vector): vector y from above.
    % U (matrix): set of extremal points of the dual problem computed by get_u(A).
    %
    % Outputs:
    % x (vector): optimal solution of the above problem.
    % val (scalar): optimal value of the above problem.
    
    tol = 1e-10;
    
    A = [A1; A2];
    n = size(A, 2);
    
    % Compute the optimal dual solution
    [val, i_max] = max(U*y);
    u_opt = U(i_max, :)';
    
    % Assign the index sets for complementarity
    I0 = abs(A'*u_opt) <= tol;
    I1 = A'*u_opt > tol;
    I2 = A'*u_opt < -tol;
    J = [true(size(A1,1), 1); u_opt(size(A1,1)+1:end) < -tol];
    
    % Use the complementarity conditions to compute the primal solution
    x = zeros(n,1);
    x(I1) = val;
    x(I2) = -val;
    x(I0) = A(J,I0) \ (y(J) - A(J,I1|I2)*x(I1|I2));
end

