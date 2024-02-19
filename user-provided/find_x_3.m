function x = find_x_3(~, d, I0, max_inf_norm)
    if isequal(I0, logical([0;1;0;1;0;1]))
        direction = [1; 1; -2];
        D_pse = [
            0.750000000000000   0.288675134594813;
            -0.750000000000000   0.288675134594813;
            0   0.288675134594813;
            ];
        idx = [1; 2];
        x0 = D_pse * d(idx);
        [s_min, s_max] = solutions_prescribed_min_norm(x0, direction, max_inf_norm);
        x = x0 + 0.5*(s_min+s_max)*direction;
    elseif isequal(I0, logical([1;0;0;1;1;0]))
        direction = [1; -1; -2];
        D_pse = [
            0.625000000000000   0.505181485540922;
            0.125000000000000   0.793856620135735;
            0.250000000000000  -0.144337567297406;
            ];
        idx = [1; 2];
        x0 = D_pse * d(idx);
        [s_min, s_max] = solutions_prescribed_min_norm(x0, direction, max_inf_norm);
        x = x0 + 0.5*(s_min+s_max)*direction;
    elseif isequal(I0, logical([1;1;1;0;0;0]))
        direction = [1; 1; 2];
        D_pse = [
            0.625000000000000  -0.505181485540922;
            -0.125000000000000   0.793856620135735;
            -0.250000000000000  -0.144337567297406;
            ];
        idx = [1; 2];
        x0 = D_pse * d(idx);
        [s_min, s_max] = solutions_prescribed_min_norm(x0, direction, max_inf_norm);
        x = x0 + 0.5*(s_min+s_max)*direction;
    elseif isequal(I0, logical([0;0;1;0;1;1]))
        direction = [1; 1; 1];
        D_pse = [
            -0.500000000000000  -0.288675134594813;
            0.500000000000000  -0.288675134594813;
            0   0.577350269189626;
            ];
        idx = [1; 2];
        x0 = D_pse * d(idx);
        [s_min, s_max] = solutions_prescribed_min_norm(x0, direction, max_inf_norm);
        x = x0 + 0.5*(s_min+s_max)*direction;
    elseif length(I0) == 6
        throw('The case above is not handled.')
    else
        x = [];
    end
end

