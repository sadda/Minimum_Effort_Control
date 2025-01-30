classdef Pars < handle
    properties
        A;
        n;
        A_case;
        A_subcase = {};
        zero_columns;
        multiples;
        multiples_counts;
        U;
        D = {};
        a = {};
        I = {};
        J = {};
        K = {};
        idx = {};
        d_vec = {};
        pars_subset = {};
    end

    methods
        function self = Pars(A)
            tol = 1e-10;
            self.A = A;
            self.n = size(A, 2);

            % Find zeros columns
            zero_columns = vecnorm(A, 2, 1) <= tol;
            A = A(:, ~zero_columns);

            % Find columns which are multiples of each other
            [multiples, multiples_counts] = find_column_multiples(A, tol);
            % Multiply the columns which are multiplied
            for i = 1:size(multiples_counts,1)
                k = multiples_counts(i, 1);
                A(:,k) = A(:,k) * multiples_counts(i, 2);
            end
            % Remove columns which are multiples
            A(:, multiples(:,2)) = [];

            self.zero_columns = zero_columns;
            self.multiples = multiples;
            self.multiples_counts = multiples_counts;

            % TODO: implement zero columns and multiples
            [m, n] = size(A);
            rank_A = rank(A);
            if rank_A ~= m
                error('Matrix does not have linearly independent rows');
            end
            if m == 1 % Theoretically this should never happen in multiples work correctly
                self.A_case = 1;
                self.D = 1/sum(abs(A))*ones(n, 1).*sign(A');
            elseif n == m
                self.A_case = 1;
                self.D = inv(A);
            elseif n == m + 1
                self.A_case = 2;
                self.D = A' / (A*A');
                self.a = null(A);
                if norm(self.a - mean(self.a)) < tol
                    self.A_subcase = 1;
                else
                    self.A_subcase = 2;
                end
            else
                % TODO: remove symmetric elements
                self.A_case = 3;
                U = get_u(A);
                for i = 1:size(U, 1)
                    ATu = A'*U(i,:)';
                    self.I{i} = (ATu <= tol) & (ATu >= -tol);
                    self.J{i} = ATu > tol;
                    self.K{i} = ATu < -tol;

                    A_I = A(:, self.I{i});
                    [A_I, self.idx{i}] = linearly_independent_rows(A_I);
                    [m_I, n_I] = size(A_I);
                    if rank(A_I) ~= m_I
                        error('Something is wrong. Matrix should have had linearly independent rows');
                    end

                    if n_I == m_I
                        self.A_subcase{i} = 1;
                        self.D{i} = inv(A_I);
                    elseif n_I == m_I + 1
                        self.A_subcase{i} = 2;
                        self.D{i} = A_I' / (A_I*A_I');
                        self.a{i} = null(A_I);
                    else
                        self.A_subcase{i} = 3;
                        self.pars_subset{i} = Pars(A_I);
                    end
                    self.d_vec{i} = sum(A(:, self.J{i}), 2) - sum(A(:, self.K{i}), 2);
                    self.U = U;
                end
            end
        end

        function x_new = expand_solution(self, x)
            % Expand the solution to the original space
            idx = find(~self.zero_columns);
            idx = idx(setdiff(1:length(idx), self.multiples(:,2)));
            x_new = zeros(self.n, 1);
            x_new(idx) = x;
            % Distribute the values into the multiples columns
            for k = 1:size(self.multiples,1)
                x_new(self.multiples(k,2)) = x_new(self.multiples(k,1)) * sign(self.multiples(k,3));
            end
        end

        function x_new = shrink_solution(self, x)
            % Shrink the solution to the reduced space
            idx = find(~self.zero_columns);
            idx = idx(setdiff(1:length(idx), self.multiples(:,2)));
            x_new = x(idx);
        end
    end
end



function [multiples, multiples_counts] = find_column_multiples(A, tol)
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