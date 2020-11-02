function [x, val, z, u] = Solve_Ours4(A, y, z, u, B, C)
    
    [~, n] = size(A);
    
    x0 = A \ y;
    
    rho = 0.5;
    rhoInv = 1/rho;
    
    for i=1:10000
        s = C * (z - x0 + u);
        
        c     = x0 + B*s - u;
        c_abs = abs(c);
        
        c_sort = sort(c_abs, 'descend');
        
        exp_sum = 0;
        z_max = Inf;
        terminated = false;
        for j=2:n
            exp_sum = exp_sum + (j-1)*(c_sort(j-1) - c_sort(j));
            if exp_sum >= rhoInv
                terminated = true;
                z_max = c_sort(j) + (exp_sum - rhoInv) / (j-1);
                break;
            end
        end
        if ~terminated
            z_max = (sum(c_sort) - rhoInv) / n;
        end
        z  = c;
        ii = c_abs >= z_max;
        z(ii) = sign(c(ii))*z_max;
        
        u_add = z - x0 - B*s;
        
        %         for j=1:n+1
        %             if j==n+1 || sum(max(c_sort - c_sort(j), 0)) >= 1/rho
        %                 break;
        %             end
        %         end
        %         z_max  = (sum(c_sort(1:j-1)) - 1/rho) / (j-1);
        %         z_old      = min(z_max,c_abs);
        %         z_old(c<0) = -z_old(c<0);
        %         norm(z - z_old)
        
        if i > 1 && norm(u_add) <= 1e-6 && norm(s - s_prev) <= 1e-6
            break;
        end
        
        u = u + u_add;
        
        s_prev = s;
        
        
    end
    
    x   = x0 + B*s;
    val = max(abs(x));
    
end


