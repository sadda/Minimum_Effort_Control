classdef Pars < handle
    properties
        A_original;
        m;
        n;
        zero_columns;
        multiples;
        multiples_counts;        
        A_case;
        A_subcase = {};
        U;
        I = {};
        J = {};
        K = {};
        D = {};
        D_idx = {};
        a = {};
        d_vec = {};
        pars_subset = {};
        tol = 1e-10;
    end

    methods
        function self = Pars(A, remove_same_columns, tol)
            if nargin < 2 || isempty(remove_same_columns)
                remove_same_columns = true;
            end
            if nargin >= 3
                self.tol = tol;
            end
            self.add_data(A, remove_same_columns);
        end

        function add_data(self, A, remove_same_columns)
            A = self.modify_input_matrix(A, remove_same_columns);

            rank_A = rank(A);
            if rank_A ~= self.m
                error('Matrix does not have linearly independent rows');
            end
            if self.m == 1
                self.A_case = 1;
                self.D = 1/sum(abs(A))*ones(self.n, 1).*sign(A');
            elseif self.n == self.m
                self.A_case = 1;
                self.D = inv(A);
            elseif self.n == self.m + 1
                self.A_case = 2;
                self.D = A' / (A*A');
                self.a = null(A);
                if norm(self.a - mean(self.a)) < self.tol
                    self.A_subcase = 1;
                else
                    self.A_subcase = 2;
                end
            else
                % TODO: remove symmetric elements
                self.A_case = 3;
                self.U = get_u(A);
                for i = 1:size(self.U, 1)
                    ATu = A'*self.U(i,:)';
                    self.I{i} = (ATu <= self.tol) & (ATu >= -self.tol);
                    self.J{i} = ATu > self.tol;
                    self.K{i} = ATu < -self.tol;

                    A_I = A(:, self.I{i});
                    [A_I, self.D_idx{i}] = linearly_independent_rows(A_I);
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
                end
            end
        end

        function A = modify_input_matrix(self, A, remove_same_columns)
            self.A_original = A;
            [self.m, self.n] = size(A);

            % Find zeros columns
            self.zero_columns = vecnorm(A, 2, 1) <= self.tol;
            A = A(:, ~self.zero_columns);

            if remove_same_columns
                % Find columns which are multiples of each other
                [self.multiples, self.multiples_counts] = find_column_multiples(A, self.tol);
                % Multiply the columns which are multiplied
                for i = 1:size(self.multiples_counts,1)
                    k = self.multiples_counts(i, 1);
                    A(:,k) = A(:,k) * self.multiples_counts(i, 2);
                end
                % Remove columns which are multiples
                A(:, self.multiples(:,2)) = [];
            else
                self.multiples = zeros(0, 3);
                self.multiples_counts = zeros(0, 2);
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

function [X_idx, idx] = linearly_independent_rows(X, tol)
    if ~nnz(X)
        error('Zero matrix');
    end
    if nargin < 2
        tol=1e-10;
    end

    % QR decomposition (transposition for passing from columns to rows)
    [~, R, E] = qr(X', 0);
    if ~isvector(R)
        diagr = abs(diag(R));
    else
        diagr = abs(R(1));
    end

    % Rank estimation
    r = find(diagr >= tol*diagr(1), 1, 'last');
    idx = sort(E(1:r));
    X_idx = X(idx,:);
end
