classdef Solver < handle
    properties
        pars
        find_x
        tol
    end

    methods
        function self = Solver(pars, find_x, tol)
            if nargin < 3
                tol = 1e-10;
            end
            if nargin < 2 || isequal(find_x, [])
                find_x = @(varargin) [];
            end

            self.pars = pars;
            self.find_x = find_x;
            self.tol = tol;
        end

        function [x, optimal_value] = min_effort(self, y)
            % min_effort Our algorithm for solving the minimum effort problem
            %    minimize     ||x||_infty
            %    subject to   self.pars.A_original * x = y
            %
            % Inputs:
            % self.pars (struct): structure of internal data obtained from get_u(A).
            % y (vector): vector y from above.
            % find_x (function, optional): upon calling find_x(D, d, I0) solves.
            %
            % Outputs:
            % x (vector): optimal solution of the above problem.
            % optimal_value (scalar): optimal value of the above problem.
            % self.pars (struct): contains information of how solutions were computed.

            [x, I0, optimal_value, D, d] = self.solve_optimal_value(y);

            % We need to solve the following equation for x(I0)
            % D * x(I0) = d;
            rank_D = rank(D);
            x_user = self.find_x(self, D, d, I0, optimal_value);
            if ~isequal(x_user, [])
                % Use user-provided solution in present
                x(I0) = x_user;
                self.pars.solution_part(I0, 'user-specified solution');
            elseif rank_D == size(D, 2)
                % The simple case when it is well-conditioned
                x(I0) = D \ d;
                self.pars.solution_part(I0, 'unique solution', 1, 'D', D);
            else
                % Reduce the system to the necessary equations only. Rank of D is still m.
                [~, idx] = rref(D');
                D = D(idx, :);
                d = d(idx);
                if size(D,1)+1 == size(D,2)
                    % Solve n*(n+1) system
                    [x0, direction, s_min, s_max] = self.solve_n_n_plus_one_all_solutions(D, d, optimal_value);
                    x(I0) = x0 + 0.5*(s_min+s_max)*direction;
                    self.pars.solution_part(I0, 'n*(n+1) system solution', 4, 'idx', idx, 'D', D, 'D_pse', D'/(D*D'), 'direction', direction, 'x0', x0, 's_min', s_min, 's_max', s_max);
                else
                    % Try l2 solution with reduced ranks
                    x(I0) = D' * ((D * D') \ d);
                    if max(abs(x)) <= optimal_value + self.tol
                        self.pars.solution_part(I0, 'l2 solution', 1, 'D', D);
                    else
                        % TODO: does not work for find_x
                        pars_subsystem = Pars(D);
                        solver_subsystem = Solver(pars_subsystem);
                        x(I0) = solver_subsystem.min_effort(d);
                        self.pars.solution_part(I0, 'get_u solution', 1, 'D', D);
                    end
                end
            end

            x = self.expand_solution(x);
            self.check_solution_quality(x, y, optimal_value);
        end

        function [x, optimal_value] = min_effort_user_provided(self, y)
            [x, I0, optimal_value, D, d] = self.solve_optimal_value(y);
            x(I0) = self.find_x(self, D, d, I0, optimal_value);
            x = self.expand_solution(x);
            self.check_solution_quality(x, y, optimal_value);
        end

        function [x, I0, optimal_value, D, d] = solve_optimal_value(self, y)
            % Compute the optimal dual solution
            [optimal_value, i_max] = max(self.pars.U*y);
            u_opt = self.pars.U(i_max, :)';

            % Assign the index sets for complementarity
            A = self.pars.A;
            [~, n] = size(A);
            I0 = abs(A'*u_opt) <= self.tol;
            I1 = A'*u_opt > self.tol;
            I2 = A'*u_opt < -self.tol;

            % Use the complementarity conditions to compute the primal solution
            x = zeros(n,1);
            x(I1) = optimal_value;
            x(I2) = -optimal_value;

            D = A(:,I0);
            d = y - A(:,I1|I2)*x(I1|I2);
        end

        function x = solve_n_n_plus_one(self, A, b)
            [m, n] = size(A);
            rank_A = rank(A);
            assert(m + 1 == n, "Matrix A must have shape (m, m+1).");
            assert(rank_A == m, "Matrix A must have rank m.");

            % Find the kernel and a particular solution
            d = null(A);
            x0 = A \ b;

            % Find the solution
            x = self.solution_min_norm(x0, d);
        end

        function x_opt = solution_min_norm(self, x0, d)
            % Finds x = x0+s*d with minimal l_infty norm.
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

        function [x0, d, s_min, s_max] = solve_n_n_plus_one_all_solutions(self, A, b, max_inf_norm)
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
            [s_min, s_max] = self.solutions_prescribed_min_norm(x0, d, max_inf_norm);
        end

        function [s_min, s_max] = solutions_prescribed_min_norm(self, x0, d, max_inf_norm)
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

        function x_expanded = expand_solution(self, x)
            % Expand the solution to the original space
            multiples = self.pars.multiples;
            idx = find(~self.pars.zero_columns);
            idx = idx(setdiff(1:length(idx), multiples(:,2)));
            x_expanded = zeros(self.pars.n,1);
            x_expanded(idx) = x;
            % Distribute the values into the multiples columns
            for k = 1:size(multiples,1)
                x_expanded(multiples(k,2)) = x_expanded(multiples(k,1)) * sign(multiples(k,3));
            end
        end

        function check_solution_quality(self, x, y, optimal_value)
            % Check for solution optimality
            if norm(self.pars.A_original*x-y) > self.tol
                throw("Problem was not solved");
            end
            if max(abs(x)) > optimal_value + self.tol
                warning("COMPUTATIONS FAILED. SOLUTION IS SOBOPTIMAL.");
            end
        end
    end
end

