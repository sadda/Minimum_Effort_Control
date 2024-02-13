function x = solve_n_n_plus_one(A, b)
    [m, n] = size(A);
    rank_A = rank(A);
    assert(m + 1 == n, "Matrix A must have shape (m, m+1).");
    assert(rank_A == m, "Matrix A must have rank m.");

    % Find the kernel and a particular solution
    d = null(A);
    x0 = A \ b;

    % Find the solution
    x = solution_min_norm(x0, d);
end