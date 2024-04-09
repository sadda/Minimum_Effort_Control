function x = find_x_3(solver, D, d, I0, max_inf_norm)
    if isequal(I0, logical([0;1;0;1;0;1]))
        direction = [1; 1; -2];
        D_pse = [
            0.750000000000000   0.288675134594813;
            -0.750000000000000   0.288675134594813;
            0   0.288675134594813;
            ];
        idx = [1; 2];
        x = find_solution(solver, I0, D, D_pse, d, idx, direction, max_inf_norm);
    elseif isequal(I0, logical([1;0;0;1;1;0]))
        direction = [1; -1; -2];
        D_pse = [
            0.625000000000000   0.505181485540922;
            0.125000000000000   0.793856620135735;
            0.250000000000000  -0.144337567297406;
            ];
        idx = [1; 2];
        x = find_solution(solver, I0, D, D_pse, d, idx, direction, max_inf_norm);
    elseif isequal(I0, logical([1;1;1;0;0;0]))
        direction = [1; 1; 2];
        D_pse = [
            0.625000000000000  -0.505181485540922;
            -0.125000000000000   0.793856620135735;
            -0.250000000000000  -0.144337567297406;
            ];
        idx = [1; 2];
        x = find_solution(solver, I0, D, D_pse, d, idx, direction, max_inf_norm);
    elseif isequal(I0, logical([0;0;1;0;1;1]))
        direction = [1; 1; 1];
        D_pse = [
            -0.500000000000000  -0.288675134594813;
            0.500000000000000  -0.288675134594813;
            0   0.577350269189626;
            ];
        idx = [1; 2];
        x = find_solution(solver, I0, D, D_pse, d, idx, direction, max_inf_norm);
    elseif isequal(I0, logical([1;0;0;0;0;1]))
        D_pse = [
            0.750000000000000   0;
            0   0.866025403784439;
            ];
        idx = [1; 2];
        x = D_pse * d(idx);
        add_counter_unique(solver, I0, D, D_pse, idx);
    elseif isequal(I0, logical([0;1;0;0;1;0]))
        D_pse = [
            0.375000000000000   0.649519052838329
            0.750000000000000  -0.433012701892219
            ];
        idx = [1; 2];
        x = D_pse * d(idx);
        add_counter_unique(solver, I0, D, D_pse, idx);
    elseif isequal(I0, logical([0;0;1;1;0;0]))
        D_pse = [
            -0.750000000000000  -0.433012701892219
            -0.375000000000000   0.649519052838329
            ];
        idx = [1; 2];
        x = D_pse * d(idx);
        add_counter_unique(solver, I0, D, D_pse, idx);
    elseif length(I0) == 6
        throw('The case above is not handled.')
    else
        x = [];
    end
end

function x = find_solution(solver, I0, D, D_pse, d, idx, direction, max_inf_norm)
    x0 = D_pse * d(idx);
    [s_min, s_max] = solver.n_n_plus_one_all_solutions(x0, direction, max_inf_norm);
    s = 0.5*(s_min+s_max);
    x = x0 + s*direction;
    add_counter_other(solver, I0, D, D_pse, idx, direction, x0, s, s_min, s_max)
end

function add_counter_unique(solver, I0, D, D_pse, idx)
    solver.pars.increase_counter(I0, 'unique solution', 3, 'D', D, 'D_pse', D_pse, 'D_idx', idx);
end

function add_counter_other(solver, I0, D, D_pse, idx, v, x0, s, s_min, s_max)
    solver.pars.increase_counter(I0, 'n*(n+1) system solution', 4, 'D', D, 'D_pse', D_pse, 'D_idx', idx, 'D_v', v, 'x0', x0, 's', s, 's_min', s_min, 's_max', s_max);
end
