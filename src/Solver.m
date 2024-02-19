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
            % y (vector): vector y from above.
            %
            % Outputs:
            % x (vector): optimal solution of the above problem.
            % optimal_value (scalar): optimal value of the above problem.

            [x, I0, optimal_value, D, d] = self.duality_solution(y);

            % We need to solve the following equation for x(I0)
            % D * x(I0) = d;
            x_user = self.find_x(self, D, d, I0, optimal_value);
            if ~isequal(x_user, [])
                % Use user-provided solution in present
                x(I0) = x_user;
                self.pars.increase_counter(I0, 'user-specified solution');
            elseif rank(D) == size(D, 2)
                % The simple case when it is well-conditioned
                x(I0) = D \ d;
                [~, idx] = rref(D');
                self.pars.increase_counter(I0, 'unique solution', 3, 'D', D, 'D_pse', inv(D(idx,:)), 'D_idx', idx);
            else
                % Reduce the system to the necessary equations only. Rank of D is still m.
                [~, idx] = rref(D');
                D = D(idx, :);
                d = d(idx);
                if size(D,1)+1 == size(D,2)
                    % Solve n*(n+1) system
                    [x0, v, s_min, s_max] = self.n_n_plus_one_all_solutions_matrix_form(D, d, optimal_value);
                    x(I0) = x0 + 0.5*(s_min+s_max)*v;
                    self.pars.increase_counter(I0, 'n*(n+1) system solution', 4, 'D', D, 'D_pse', D'/(D*D'), 'D_idx', idx, 'D_v', v, 'x0', x0, 's_min', s_min, 's_max', s_max);
                else
                    % Try l2 solution with reduced ranks
                    x(I0) = D' * ((D * D') \ d);
                    if max(abs(x)) <= optimal_value + self.tol
                        self.pars.increase_counter(I0, 'l2 solution', 1, 'D', D);
                    else
                        % TODO: does not work for find_x
                        pars_subsystem = Pars(D);
                        solver_subsystem = Solver(pars_subsystem);
                        x(I0) = solver_subsystem.min_effort(d);
                        self.pars.increase_counter(I0, 'get_u solution', 1, 'D', D);
                    end
                end
            end

            x = self.expand_solution(x);
            self.check_solution_quality(x, y, optimal_value);
        end

        function [x, optimal_value] = min_effort_user_provided(self, y)
            [x, I0, optimal_value, D, d] = self.duality_solution(y);
            x(I0) = self.find_x(self, D, d, I0, optimal_value);
            x = self.expand_solution(x);
            self.check_solution_quality(x, y, optimal_value);
        end

        function [x, I0, optimal_value, D, d] = duality_solution(self, y)
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

        function x = n_n_plus_one_min_norm_solution_matrix_form(self, A, b)
            self.check_n_n_plus_one(A)

            % Find the kernel and a particular solution
            d = null(A);
            x0 = A \ b;

            % Find the solution
            x = self.n_n_plus_one_min_norm_solution(x0, d);
        end

        function x_opt = n_n_plus_one_min_norm_solution(~, x0, d)
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

        function [x0, v, s_min, s_max] = n_n_plus_one_all_solutions_matrix_form(self, A, b, max_inf_norm)
            self.check_n_n_plus_one(A)

            % Find the kernel and a particular solution
            v = null(A);
            if v(1) ~= 0
                v = v ./ v(1);
                v = v ./ sign(v(1));
            end
            x0 = A \ b;

            % Find the solution
            [s_min, s_max] = self.n_n_plus_one_all_solutions(x0, v, max_inf_norm);
        end

        function [s_min, s_max] = n_n_plus_one_all_solutions(~, x0, v, max_inf_norm)
            % Finds all solution x = x0+s*v with l_infty norm below max_inf_norm.
            % In the case of v=null(A) and x0=A\b when A has size n*(n+1),
            % the function finds the solution of A*x=b.

            s_min = -inf;
            s_max = inf;
            for i=1:length(v)
                if v(i) > 0
                    s_min = max(s_min, (-max_inf_norm - x0(i)) / v(i));
                    s_max = min(s_max, (max_inf_norm - x0(i)) / v(i));
                else
                    s_min = max(s_min, (max_inf_norm - x0(i)) / v(i));
                    s_max = min(s_max, (-max_inf_norm - x0(i)) / v(i));
                end
            end
        end

        function x_expanded = expand_solution(self, x)
            % Expand the solution to the original space
            multiples = self.pars.multiples;
            idx = find(~self.pars.zero_columns);
            idx = idx(setdiff(1:length(idx), multiples(:,2)));
            x_expanded = zeros(size(self.pars.A_original,2),1);
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

        function check_n_n_plus_one(~, A)
            [m, n] = size(A);
            assert(m + 1 == n, "Matrix A must have shape (m, m+1).");
            assert(rank(A) == m, "Matrix A must have rank m.");
        end
    end
end

