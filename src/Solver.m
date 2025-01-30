classdef Solver < handle
    properties
        pars
        tol = 1e-10;
    end

    methods
        function self = Solver(pars, tol)
            if nargin >= 2
                self.tol = tol;
            end
            self.pars = pars;
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

            pars = self.pars; %#ok<*PROPLC>
            if pars.A_case == 1
                x = pars.D * y;
                optimal_value = max(abs(x));
            elseif pars.A_case == 2
                if pars.A_subcase == 1
                    x = pars.D * y;
                    x = x - 0.5*(min(x)+max(x));
                    optimal_value = max(abs(x));
                elseif pars.A_subcase == 2
                    % TODO: implement
                    error('not implemented yet');
                else
                    error('case not known');
                end
            elseif pars.A_case == 3
                [optimal_value, i_max] = max(self.pars.U*y);

                I = pars.I{i_max};
                J = pars.J{i_max};
                K = pars.K{i_max};

                % Use the complementarity conditions to compute the primal solution
                x = zeros(length(I), 1);
                x(J) = optimal_value;
                x(K) = -optimal_value;

                d = y - pars.d_vec{i_max}*optimal_value;
                d = d(pars.idx{i_max});

                if pars.A_subcase{i_max} == 1
                    x(I) = pars.D{i_max}*d;
                elseif pars.A_subcase{i_max} == 2
                    % TODO: handle any a_i=0
                    u = pars.D{i_max}*d;
                    a = pars.a{i_max};
                    s_max = max(-optimal_value./abs(a) - u./a);
                    s_min = min(optimal_value./abs(a) - u./a);
                    if s_min < s_max
                        error('s_min larger than s_max');
                    end
                    x(I) = u + 0.5*(s_min+s_max)*a;
                elseif pars.A_subcase{i_max} == 3
                    solver_new = Solver(pars.pars_subset{i_max}, self.tol);
                    x(I) = solver_new.min_effort(d);
                else
                    error('case not known');
                end
            else
                error('case not known');
            end
            x = self.pars.expand_solution(x);
            self.check_solution_quality(x, y, optimal_value);
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

        function check_solution_quality(self, x, y, optimal_value)
            % Check for solution optimality
            if norm(self.pars.A*x-y) > self.tol
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

