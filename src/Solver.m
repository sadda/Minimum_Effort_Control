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

            % Specify find_x if not provided
            A = self.pars.A;
            U = self.pars.U;
            [~, n] = size(A);

            % Compute the optimal dual solution
            [optimal_value, i_max] = max(U*y);
            u_opt = U(i_max, :)';

            % Assign the index sets for complementarity
            I0 = abs(A'*u_opt) <= self.tol;
            I1 = A'*u_opt > self.tol;
            I2 = A'*u_opt < -self.tol;

            % Use the complementarity conditions to compute the primal solution
            x = zeros(n,1);
            x(I1) = optimal_value;
            x(I2) = -optimal_value;

            D = A(:,I0);
            d = y - A(:,I1|I2)*x(I1|I2);

            % We need to solve the following equation for x(I0)
            % D * x(I0) = d;
            rank_D = rank(D);
            x_user = self.find_x(D, d, I0, optimal_value);
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
                    [x0, direction, s_min, s_max] = solve_n_n_plus_one_all_solutions(D, d, optimal_value);
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

            % Expand the solution to the original space
            multiples = self.pars.multiples;
            x2 = zeros(self.pars.n,1);
            idx = find(~self.pars.zero_columns);
            idx = idx(setdiff(1:length(idx), multiples(:,2)));
            x2(idx) = x;
            % Distribute the values into the multiples columns
            for k = 1:size(multiples,1)
                x2(multiples(k,2)) = x2(multiples(k,1)) * sign(multiples(k,3));
            end
            x = x2;

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

