function x = find_x_3(D, d, I0)
    tol = 1e-10;
    if isequal(I0, logical([0;0;1;0;0;1;0;0;1]))
        D = D(1:2, :);
        d = d(1:2);

        D_pse = D'*inv(D*D');
        x = D_pse * d;
        x = x - 0.5*(min(x)+max(x));
    elseif isequal(I0, logical([0;1;0;1;0;0;1;1;1]))
        D(:, [3, 4]) = 2*D(:, [3, 4]);
        D = D(1:2, [3, 4, 5]);
        d = d(1:2);
        D_pse = D'*inv(D*D');
        x = D_pse * d;
        x1 = optimal_x(x, [1;1;2]);
        x = [-x1(1); -x1(2); x1(1); x1(2); x1(3)];
    elseif isequal(I0, logical([1;1;1;0;1;0;1;0;0]))
        D(:, [1, 2]) = 2*D(:, [1, 2]);
        D = D(1:2, [1, 2, 3]);
        d = d(1:2);
        D_pse = D'*inv(D*D');
        x = D_pse * d;
        x1 = optimal_x(x, [1;1;2]);
        x = [x1(1); x1(2); x1(3); -x1(1); -x1(2)];
    elseif isequal(I0, logical([1;0;0;1;1;1;0;1;0]))
        % Goes into the l2 solution which is fine here
        x = [];
    elseif isequal(I0, logical([0;1;0;0;0;1;1;0;0]))
        % Goes into the l2 solution which is fine here
        x = [];
    elseif length(I0) == 9
        throw('The case above is not handled.')
    else
        x = [];
    end
end

function x_opt = optimal_x(x0, d)
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