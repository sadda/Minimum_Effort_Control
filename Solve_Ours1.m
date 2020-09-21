function [x, val] = Solve_Ours1(A, y)
    
    mn = size(A);
    m  = mn(1);
    n  = mn(2);
    
    %% Primal problem (no linprog)
    
    % The solution set of Ax=y is equivalent to finding a particular solution
    % x0 and considering x0 + ker(A), where ker(A) is the kernel of a matrix.
    
    % If there is one degree of freedom, this corresponds to x0+s*a.
    % Then we look for s such that either
    % x0(i)+s*a(i) = x0(j)+s*a(j) or
    % x0(i)+s*a(i) = -(x0(j)+s*a(j))
    
    
    if m ~= n-1
        % error("Works only for matrices (m,m+1).");
        x   = nan(n, 1);
        val = nan(1);
        return;
    end
    
    a  = null(A); % [1;1;1] with (A*a=0)
    x0 = A \ y;
    
    f_min = Inf;
    for i = 1:n
        for j = i+1:n
            if a(i) ~= a(j)
                s = (-x0(i)+x0(j)) / (a(i)-a(j));
                if max(abs(x0 + s*a)) <= f_min
                    f_min = max(abs(x0 + s*a));
                    s_opt = s;
                end
            end
            if a(i) ~= -a(j)
                s = (-x0(i)-x0(j)) / (a(i)+a(j));
                if max(abs(x0 + s*a)) <= f_min
                    f_min = max(abs(x0 + s*a));
                    s_opt = s;
                end
            end
        end
    end
    x   = x0 + s_opt*a;
    val = max(abs(x));
    
end

