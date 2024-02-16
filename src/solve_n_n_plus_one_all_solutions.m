function [x0, d, s_min, s_max] = solve_n_n_plus_one_all_solutions(A, b, max_inf_norm)
    [m, n] = size(A);
    rank_A = rank(A);
    assert(m + 1 == n, "Matrix A must have shape (m, m+1).");
    assert(rank_A == m, "Matrix A must have rank m.");

    % Find the kernel and a particular solution
    d = null(A);
    if d(1) ~= 0
        d = d ./ d(1);
        d = d ./ sign(d(1));
    end
    x0 = A \ b;

    % Find the solution
    [s_min, s_max] = solutions_prescribed_min_norm(x0, d, max_inf_norm);
end