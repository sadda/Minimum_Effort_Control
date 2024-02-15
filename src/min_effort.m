function [x, optimal_value, pars] = min_effort(pars, y, find_x)
    % min_effort Our algorithm for solving the minimum effort problem
    %    minimize     ||x||_infty
    %    subject to   A * x = y1
    %
    % Inputs:
    % A (matrix): matrix A from above.
    % y (vector): vector y from above.
    % U (matrix): set of extremal points of the dual problem computed by get_u(A, B).
    % find_x (function, optional): upon calling find_x(D, d, I0) solves.
    % s (struct, optional): used for monitoring where the solution is computed.
    %       the minimum effort problem D * x_out = d. Even though
    %       this procedure should handle any case, it may used to handle
    %       problematic cases in a fast way.
    %
    % Outputs:
    % x (vector): optimal solution of the above problem.
    % optimal_value (scalar): optimal value of the above problem.

    % Specify find_x if not provided
    if nargin < 3 || isequal(find_x, [])
        find_x = @(varargin) [];
    end

    A = pars.A;
    U = pars.U;
    tol = 1e-10;
    [~, n] = size(A);

    % Compute the optimal dual solution
    [optimal_value, i_max] = max(U*y);
    u_opt = U(i_max, :)';

    % Assign the index sets for complementarity
    I0 = abs(A'*u_opt) <= tol;
    I1 = A'*u_opt > tol;
    I2 = A'*u_opt < -tol;

    % Use the complementarity conditions to compute the primal solution
    x = zeros(n,1);
    x(I1) = optimal_value;
    x(I2) = -optimal_value;

    D = A(:,I0);
    d = y - A(:,I1|I2)*x(I1|I2);

    % We need to solve the following equation for x(I0)
    % D * x(I0) = d;
    rank_D = rank(D);
    x_user = find_x(D, d, I0);
    if ~isequal(x_user, [])
        % Use user-provided solution in present
        x(I0) = x_user;
        pars = solution_part(pars, I0, 'user-specified');
    elseif rank_D == n
        % The simple case when it is well-conditioned
        x(I0) = D \ d;
        pars = solution_part(pars, I0, 'well-conditioned');
    else
        % Remove zero columns
        D(:, sum(abs(D),1) == 0) = [];
        % Reduce the system to the necessary equations only. Rank of D is still m.
        [~, idx] = rref(D');
        %D = D(1:rank_D,:);
        %d = d(1:rank_D);
        D = D(idx, :);
        d = d(idx);
        [m, n] = size(D);
        % Try l2 solution with reduced ranks
        x(I0) = D' * ((D * D') \ d);
        if max(abs(x)) > optimal_value + tol
            if m+1 == n
                % Solve n*(n+1) system
                x(I0) = solve_n_n_plus_one(D, d);
                pars = solution_part(pars, I0, 'n*(n+1) system', 'row_indices', idx);
            else
                U_D = get_u(D);
                x(I0) = min_effort(D, d, U_D);
                pars = solution_part(pars, I0, 'get_u');
            end
        else
            pars = solution_part(pars, I0, 'l2 with reduced', 'row_indices', idx);
        end
    end

    % Expand the solution to the original space
    multiples = pars.multiples;
    x2 = zeros(pars.n,1);
    idx = find(~pars.zero_columns);
    idx = idx(setdiff(1:length(idx), multiples(:,2)));
    x2(idx) = x;
    % Distribute the values into the multiples columns
    for k = 1:size(multiples,1)
        x2(multiples(k,2)) = x2(multiples(k,1)) * sign(multiples(k,3));
    end
    x = x2;

    % Check for solution optimality
    if norm(pars.A_original*x-y) > tol
        throw("Problem was not solved");
    end
    if max(abs(x)) > optimal_value + tol
        warning("COMPUTATIONS FAILED. SOLUTION IS SOBOPTIMAL.");
    end
end



function pars = solution_part(pars, I0, text, varargin)
    if ~isfield(pars, "analysis")
        pars.analysis = struct('I0', I0, 'text', text, 'count', 1);
        for i = 1:length(varargin)/2
            pars.analysis = {create_new_struct(I0, text, varargin{:})};
        end
    else
        found = false;
        for i = 1:length(pars.analysis)
            if isequal(pars.analysis{i}.I0, I0) && isequal(pars.analysis{i}.text, text)
                found = true;
                break
            end
        end
        if found
            pars.analysis{i}.count = pars.analysis{i}.count + 1;
        else
            pars.analysis = [pars.analysis; create_new_struct(I0, text, varargin{:})];
        end
    end
end

function output = create_new_struct(I0, text, varargin)
    output = struct('I0', I0, 'text', text, 'count', 1);
    for i = 1:length(varargin)/2
        output.(varargin{2*i-1}) = varargin{2*i};
    end
end
