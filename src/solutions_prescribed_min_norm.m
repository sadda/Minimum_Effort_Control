function [s_min, s_max] = solutions_prescribed_min_norm(x0, d, max_inf_norm)
    % Finds all solution x = x0+s*d with l_infty norm below max_inf_norm.
    n = length(d);
    s_min = -inf;
    s_max = inf;
    for i=1:n
        if d(i) > 0
            s_min = max(s_min, (-max_inf_norm - x0(i)) / d(i));
            s_max = min(s_max, (max_inf_norm - x0(i)) / d(i));
        else
            s_min = max(s_min, (max_inf_norm - x0(i)) / d(i));
            s_max = min(s_max, (-max_inf_norm - x0(i)) / d(i));
        end
    end
end