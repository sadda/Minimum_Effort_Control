function [x, val] = Solve_Ours2(A, y)
    
    mn = size(A);
    m  = mn(1);
    n  = mn(2);
    
    a  = ones(size(A,2),1);
    if norm(A*a) >= 1e-10 || m ~= n-1
        % error('Does not work in this setting');
        x   = nan(n, 1);
        val = nan(1);
        return;
    end
    
    x0    = A \ y;
    
    s_opt = -0.5*(max(x0) + min(x0));
    x     = x0 + s_opt;
    val   = max(abs(x));
    
end