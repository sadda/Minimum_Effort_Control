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

        function [x, optimal_value] = l2_solution(self, y)
            % l2_solution Our algorithm for solving the l2 solution
            %    minimize     ||x||_2
            %    subject to   self.pars.A_original * x = y
            %
            % Inputs:
            % y (vector): vector y from above.
            %
            % Outputs:
            % x (vector): optimal solution of the above problem.
            % optimal_value (scalar): optimal value of the above problem.

            if self.pars.m_B > 0
                error('l2 solution only works for B = [].');
            end
            x = self.pars.D_l2 * y;
            optimal_value = max(abs(x));
        end

        function [x, optimal_value] = min_effort(self, y)
            % min_effort Our algorithm for solving the minimum effort problem
            %    minimize     ||x||_infty
            %    subject to   self.pars.A_original * x = y
            %                 self.pars.B_original * x <= y
            %
            % Inputs:
            % y (vector): vector y from above.
            %
            % Outputs:
            % x (vector): optimal solution of the above problem.
            % optimal_value (scalar): optimal value of the above problem.

            if self.pars.A_case == 1
                x = self.pars.D * y;
                optimal_value = max(abs(x));
            elseif self.pars.A_case == 2
                if self.pars.A_subcase == 1
                    x = self.pars.D * y;
                    x = x - 0.5*(min(x)+max(x));
                elseif self.pars.A_subcase == 2
                    x = self.solve_n_n_plus_one(self.pars.D * y, self.pars.a);
                end
                optimal_value = max(abs(x));
            elseif self.pars.A_case == 3
                % Find index of the optimal dual solution
                [optimal_value, i_max] = max(self.pars.U*y);
                I = self.pars.I{i_max};
                J = self.pars.J{i_max};
                K = self.pars.K{i_max};

                % Prescibe primal solution on the non-active indices J and K
                x = zeros(length(I), 1);
                x(J) = optimal_value;
                x(K) = -optimal_value;

                if self.pars.A_subcase{i_max} ~= 0            
                    % Modify the right-hand side for solving A_I*x(I) = d
                    d = y - self.pars.d_vec{i_max}*optimal_value;
                    d_A = d(1:self.pars.m_A);
                    d_B = d(self.pars.m_A+1:end);
                    d_A = [d_A; d_B(self.pars.v_idx{i_max})];
                    d_A = d_A(self.pars.D_idx{i_max});
                    d_B = d_B(~self.pars.v_idx{i_max});
                    d = [d_A; d_B];
                end
                if self.pars.A_subcase{i_max} == 1
                    x(I) = self.pars.D{i_max}*d_A;
                elseif self.pars.A_subcase{i_max} == 2
                    % Get a and u for non-zero components of a
                    u = self.pars.D{i_max}*d;
                    a = self.pars.a{i_max};
                    idx = abs(a) > self.tol;
                    % Find s_max and s_min
                    s_max = max(-optimal_value./abs(a(idx)) - u(idx)./a(idx));
                    s_min = min(optimal_value./abs(a(idx)) - u(idx)./a(idx));
                    if s_min < s_max - self.tol
                        error('s_min larger than s_max');
                    end
                    x(I) = u + 0.5*(s_min+s_max)*a;
                elseif self.pars.A_subcase{i_max} == 3
                    solver_new = Solver(self.pars.pars_subset{i_max}, self.tol);
                    x(I) = solver_new.min_effort(d);
                end
            end
            x = self.pars.expand_solution(x);
            self.check_solution_quality(x, y, optimal_value);
        end

        function x_opt = solve_n_n_plus_one(~, x0, d)
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

        function check_solution_quality(self, x, y, optimal_value)
            % Check for solution optimality
            if norm(self.pars.A_original*x-y(1:self.pars.m_A)) > self.tol
                throw("Problem was not solved");
            end
            if norm(max(self.pars.B_original*x-y(self.pars.m_A+1:end), 0)) > self.tol
                throw("Problem was not solved");
            end
            if max(abs(x)) > optimal_value + self.tol
                warning("COMPUTATIONS FAILED. SOLUTION IS SOBOPTIMAL.");
            end
        end
    end
end

