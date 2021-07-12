function [x, val] = min_effort_linprog(A, y)
    
    % mnmz ||x||_infty
    % s.t. Ax = y
    
    % Primal problem (equivalent)
    % mnmz z
    % s.t. Ax = y
    %      z  >= x_i
    %      z  >= -x_i
    
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

