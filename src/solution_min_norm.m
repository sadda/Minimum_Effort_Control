function x_opt = solution_min_norm(x0, d)
    % Finds x = x0+s*d with minimal l_infty norm.
    n = length(d);
    val_opt = inf;
    for i=1:n
        for j=i+1:n
            if d(i) ~= d(j)
                s = (x0(j) - x0(i)) / (d(i) - d(j));
                val = max(abs(x0 + s*d));
                if val < val_opt
                    val_opt = val;
                    x_opt = x0 + s*d;
                end
            end
            if d(i) ~= -d(j)
                s = -(x0(j) + x0(i)) / (d(i) + d(j));
                val = max(abs(x0 + s*d));
                if val < val_opt
                    val_opt = val;
                    x_opt = x0 + s*d;
                end
            end
        end
    end
end