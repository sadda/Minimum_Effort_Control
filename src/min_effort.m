function [x, optimal_value, s] = min_effort(A, y, U, find_x, s)
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
    if nargin < 5
        s = [];
    end
    if nargin < 4 || isequal(find_x, [])
        find_x = @(varargin) [];
    end    

    tol = 1e-10;
    [m, n] = size(A);

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
        s = solution_part(s, I0, 'user-specified');
    elseif rank_D == n
        % The simple case when it is well-conditioned
        x(I0) = D \ d;
        s = solution_part(s, I0, 'well-conditioned');
    else
        % Remove zero columns
        D(:, sum(abs(D),1) == 0) = [];
        % Reduce the system to the necessary equations only. Rank of D is still m.
        D = D(1:rank_D,:);
        d = d(1:rank_D);
        [m, n] = size(D);
        % Try l2 solution with reduced ranks
        x(I0) = D' * ((D * D') \ d);
        if max(abs(x)) > optimal_value + tol
            if m+1 == n
                % Solve n*(n+1) system
                x(I0) = solve_n_n_plus_one(D, d);
                s = solution_part(s, I0, 'n*(n+1) system');
            else
                % Find columns which are multiples of each other
                [multiples, multiples_counts] = find_column_multiples(D);
                % Multiply the columns which are multiplied
                D_modified = D;
                for i = 1:size(multiples_counts,1)
                    k = multiples_counts(i, 1);
                    D_modified(:,k) = D_modified(:,k) * multiples_counts(i, 2);
                end
                % Remove columns which are multiples
                D_modified(:, multiples(:,2)) = [];

                if size(D_modified,1) + 1 == size(D_modified,2)
                    % Solve the n*(n+1) system
                    x0 = zeros(sum(I0),1);
                    x0(setdiff(1:length(x0), multiples(:,2))) = solve_n_n_plus_one(D_modified, d);
                    % Distribute the values into the multiples columns
                    for k = 1:size(multiples,1)
                        x0(multiples(k,2)) = x0(multiples(k,1)) * sign(multiples(k,3));
                    end
                    x(I0) = x0;
                    s = solution_part(s, I0, 'n*(n+1) system modified');
                else
                    U_D = get_u(D);
                    x(I0) = min_effort(D, d, U_D);
                    s = solution_part(s, I0, 'get_u');
                end
            end
        else
            s = solution_part(s, I0, 'l2 with reduced');
        end        
    end

    % Check for solution optimality
    if norm(A*x-y) > tol
        throw("Problem was not solved");
    end
    if max(abs(x)) > optimal_value + tol
        warning("COMPUTATIONS FAILED. SOLUTION IS SOBOPTIMAL.");
    end
end



function [multiples, multiples_counts] = find_column_multiples(A)
    tol = 1e-10;
    n = size(A, 2);
    multiples = zeros(0, 3);
    for i = 1:n
        if ~ismember(i, multiples(:,2))
            for j = i+1:n
                v1 = A(:,i);
                v2 = A(:,j);
                c = v2 \ v1;
                if norm(v1 - c*v2) <= tol
                    multiples = [multiples; [i, j, c]];
                end
            end
        end
    end
    multiples_unique = unique(multiples(:,1));
    multiples_sum = zeros(length(multiples_unique), 1);
    for i = 1:length(multiples_unique)
        multiples_sum(i) = sum(abs(multiples(multiples(:,1)==multiples_unique(i),3))) + 1;
    end
    multiples_counts = [multiples_unique, multiples_sum];
end

function s = solution_part(s, I0, text)
    if isequal(s, [])
        s = add_field(s, I0, text);
    else
        found = false;
        for i = 1:length(s)
            if isequal(s(i).I0, I0) && isequal(s(i).text, text)
                found = true;
                break
            end
        end
        if found
            s(i).count = s(i).count + 1;
        else
            s = add_field(s, I0, text);
        end
    end
end

function s = add_field(s, I0, text)
    s = [s; struct('I0', I0, 'text', text, 'count', 1)];
end

