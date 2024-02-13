function x_out = find_x_2(A_big, y, I0, J, x)
    if isequal(I0, logical([1;0;1;1;0])) && (isequal(J, logical([1;1;0;0;0;1])) || isequal(J, logical([1;1;0;1;0;0])))
        optimal_value = max(abs(x));

        t_p = A_big([1 2],I0) \ (y([1 2]) - A_big([1 2],~I0)*x(~I0));
        t_d = null(A_big(J,I0));
        if any(t_d <= 0)
            t_d = -t_d;
        end
        c_max1 = min((optimal_value-t_p) ./ t_d);
        c_min1 = max((-optimal_value-t_p) ./ t_d);
        c_max2 = (y(3) - A_big(3,~I0)*x(~I0) - A_big(3,I0)*t_p) / (A_big(3,I0)*t_d);
        c_min2 = (-y(3) - A_big(3,~I0)*x(~I0) - A_big(3,I0)*t_p) / (A_big(3,I0)*t_d);
        c_max = min(c_max1, c_max2);
        c_min = max(c_min1, c_min2);

        xs_a = x;
        xs_a(I0) = t_p + c_min*t_d;
        xs_b = x;
        xs_b(I0) = t_p + c_max*t_d;

        norm_a = norm(A_big([3 4],:)*xs_a);
        norm_b = norm(A_big([3 4],:)*xs_b);

        tol = 1e-10;
        if norm(norm_a) > norm(norm_b) + tol
            c = c_min;
        elseif norm(norm_b) > norm(norm_a) + tol
            c = c_max;
        else
            c = 0.5*(c_min+c_max);
        end
        x_out = t_p + c*t_d;
    else
        x_out = [];
    end
end