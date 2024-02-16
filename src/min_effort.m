function [x, optimal_value, pars] = min_effort(pars, y, find_x)
    % min_effort Our algorithm for solving the minimum effort problem
    %    minimize     ||x||_infty
    %    subject to   pars.A_original * x = y
    %
    % Inputs:
    % pars (struct): structure of internal data obtained from get_u(A).
    % y (vector): vector y from above.
    % find_x (function, optional): upon calling find_x(D, d, I0) solves.
    %
    % Outputs:
    % x (vector): optimal solution of the above problem.
    % optimal_value (scalar): optimal value of the above problem.
    % pars (struct): contains information of how solutions were computed.

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
        pars = solution_part(pars, I0, 'user-specified solution');
    elseif rank_D == size(D, 2)
        % The simple case when it is well-conditioned
        x(I0) = D \ d;
        pars = solution_part(pars, I0, 'unique solution', 1, 'D', D);
    else
        % Reduce the system to the necessary equations only. Rank of D is still m.
        [~, idx] = rref(D');
        D = D(idx, :);
        d = d(idx);
        if size(D,1)+1 == size(D,2)
            % Solve n*(n+1) system
            [x0, direction, s_min, s_max] = solve_n_n_plus_one_all_solutions(D, d, optimal_value);
            x(I0) = x0 + 0.5*(s_min+s_max)*direction;
            pars = solution_part(pars, I0, 'n*(n+1) system solution', 2, 'D', D, 'D_pse', D'/(D*D'), 'direction', direction, 'x0', x0, 's_min', s_min, 's_max', s_max);
        else
            % Try l2 solution with reduced ranks
            x(I0) = D' * ((D * D') \ d);
            if max(abs(x)) <= optimal_value + tol
                pars = solution_part(pars, I0, 'l2 solution', 1, 'D', D);
            else
                % TODO: does not work for find_x
                pars_subsystem = get_u(D);
                x(I0) = min_effort(pars_subsystem, d, find_x);
                pars = solution_part(pars, I0, 'get_u solution', 1, 'D', D);
            end
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
    if pars.analysis_update
        if ~isfield(pars, 'analysis')
            pars.analysis_i = pars.analysis_i + 1;
            pars.analysis = {create_new_struct(I0, text)};
            pars = update_pars(pars, 1, varargin{:});
        else
            found = false;
            for i = 1:length(pars.analysis)
                if isequal(pars.analysis{i}.I0, I0) && isequal(pars.analysis{i}.text, text)
                    found = true;
                    break
                end
            end
            if found
                pars.analysis_i = pars.analysis_i + 1;
                pars = update_pars(pars, i, varargin{:});
            else
                pars.analysis_i = pars.analysis_i + 1;
                pars.analysis = [pars.analysis; create_new_struct(I0, text)];
                pars = update_pars(pars, length(pars.analysis), varargin{:});
            end
        end
    end
end

function output = create_new_struct(I0, text)
    output = struct('I0', I0, 'text', text, 'count', 0, 'i', []);
end

function pars = update_pars(pars, i, varargin)
    pars.analysis{i}.count = pars.analysis{i}.count + 1;
    pars.analysis{i}.i = [pars.analysis{i}.i; pars.analysis_i];
    n_constant = varargin{1};
    for j = 1:n_constant
        pars.analysis{i}.(varargin{2*j}) = varargin{2*j+1};
    end
    for j = n_constant+1:(length(varargin)-1)/2
        if isfield(pars.analysis{i}, varargin{2*j})
            pars.analysis{i}.(varargin{2*j}) = [pars.analysis{i}.(varargin{2*j}); varargin{2*j+1}'];
        else
            pars.analysis{i}.(varargin{2*j}) = [varargin{2*j+1}'];
        end
    end
end

