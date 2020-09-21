function [x, val] = Solve_Ours3(A, y)
    
    mn = size(A);
    m  = mn(1);
    n  = mn(2);
    
    a  = null(A);
    aa = a'*a \ a';
    
    x0 = A \ y;
    
    % Constraint z - x0 - a*s = 0
    
    rho = 0.5;
    
    u = ones(n, 1);
    z = ones(n, 1);
    
    for i=1:100
        s = aa * (z - x0 + u);
        
        c     = - u + x0 + a*s;
        c_abs = abs(c);
        
        %         n_try  = 10001;
        %         z_try  = linspace(0,max(c_abs),n_try);
        %         f_try  = sum(max(repmat(c_abs,1,n_try) - repmat(z_try,n,1),0));
        %         [~,ii] = min(abs(f_try - 1/rho));
        %         z_max  = z_try(ii);
        %         z      = min(z_max,c_abs);
        %         z(c<0) = -z(c<0);
        
        c_sort = sort(c_abs, 'descend');
        for j=1:n+1
            if j==n+1 || sum(max(c_sort - c_sort(j), 0)) >= 1/rho
                break;
            end
        end
        z_max  = (sum(c_sort(1:j-1)) - 1/rho) / (j-1);
        z      = min(z_max,c_abs);
        z(c<0) = -z(c<0);
        
        u = u + z - x0 - a*s;
    end
    
    x   = x0 + a*s;
    val = max(abs(x));
    
end


