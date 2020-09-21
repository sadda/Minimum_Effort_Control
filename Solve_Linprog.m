function [x, val] = Solve_Linprog(A, y, primal)
    
    mn = size(A);
    m  = mn(1);
    n  = mn(2);
    
    if primal
        %% Primal problem
        
        % mnmz ||x||_infty
        % s.t. Ax = y
        
        % Primal problem (equivalent)
        % mnmz z
        % s.t. Ax = y
        %      z  >= x_i
        %      z  >= -x_i
        
        c      = zeros(n+1,1);
        c(end) = 1;
        A_ineq = [eye(n), -ones(n,1); -eye(n), -ones(n,1)];
        b_ineq = [zeros(n,1); zeros(n,1)];
        A_eq   = [A, zeros(m,1)];
        b_eq   = y;
        lb     = [];
        ub     = [];
        
        [xz, val] = linprog(c, A_ineq, b_ineq, A_eq, b_eq, lb, ub);
        x = xz(1:n);
        
    else
        %% Dual problem
        
        % mxmz y'*u
        % s.t. ||A'*u||_1 <= 1
        
        % Dual problem (equivalent)
        % mxmz y'*u
        % s.t. A'*u = zP - zM
        %      sum(zP + zM) <= 1
        %      zP, zM >= 0
        
        
        c      = zeros(m+2*n,1);
        c(1:m) = y;
        A_ineq = [zeros(1,m), ones(1,n), ones(1,n)];
        b_ineq = 1;
        A_eq   = [A', -eye(n), eye(n)];
        b_eq   = zeros(n,1);
        lb     = [-Inf(m,1); zeros(n,1); zeros(n,1)];
        ub     = [];
        
        [uz, val, ~, ~, lambda] = linprog(-c, A_ineq, b_ineq, A_eq, b_eq, lb, ub);
        
        x   = lambda.eqlin;
        val = -val;
        
        % u = uz(1:n);
    end
end

