function [x, val] = min_effort(A, y, U)
    
    % mnmz ||x||_infty
    % s.t. Ax = y
    
    % Primal problem (equivalent)
    % mnmz z
    % s.t. Ax = y
    %      z  >= x_i
    %      z  >= -x_i
    
    tol = 1e-10;
    
    n = size(A, 2);
    
    [val, i_max] = max(U*y);
    u_opt = U(i_max, :)';
    
    I0 = abs(A'*u_opt) <= tol;
    I1 = A'*u_opt > tol;
    I2 = A'*u_opt < -tol;
    
    x = zeros(n,1);
    x(I1) = val;
    x(I2) = -val;
    x(I0) = A(:,I0) \ (y - A(:,I1|I2)*x(I1|I2));
end

