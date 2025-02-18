classdef Pars < handle
    properties
        A_original;
        m_A;
        B_original;
        m_B;
        n;
        zero_columns;
        multiples;
        multiples_counts;        
        A_case;
        A_subcase = {};
        D_l2;
        U;
        I = {};
        J = {};
        K = {};
        D = {};
        D_idx = {};
        a = {};
        d_vec = {};
        v_idx = {};
        pars_subset = {};
        tol = 1e-10;
    end

    methods
        function self = Pars(A, B, remove_same_columns, tol)
            if nargin < 2 || isempty(B)
                B = zeros(0, size(A,2));
            end
            if nargin < 3 || isempty(remove_same_columns)
                remove_same_columns = true;
            end
            if nargin >= 4
                self.tol = tol;
            end
            self.add_data(A, B, remove_same_columns);
        end

        function add_data(self, A, B, remove_same_columns)
            [A, B] = self.modify_input_matrix(A, B, remove_same_columns);
            self.check_A()

            self.D_l2 = A' / (A*A');
            if self.m_A == 1 && self.m_B == 0
                self.A_case = 1;
                self.D = 1/sum(abs(A))*ones(self.n, 1).*sign(A');
            elseif self.n == self.m_A
                self.A_case = 1;
                self.D = inv(A);
            elseif self.n == self.m_A + 1 && self.m_B == 0
                % TODO: implement for B
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
                self.U = get_u(A, B);
                for i = 1:size(self.U, 1)
                    u = self.U(i, 1:self.m_A)';
                    v = self.U(i, self.m_A+1:end)';
                    ABTu = A'*u + B'*v;
                    self.I{i} = (ABTu <= self.tol) & (ABTu >= -self.tol);
                    self.J{i} = ABTu > self.tol;
                    self.K{i} = ABTu < -self.tol;

                    if sum(self.I{i}) == 0
                        self.A_subcase{i} = 0;
                        continue;
                    end

                    self.v_idx{i} = v < -self.tol;
                    A_I = A(:, self.I{i});
                    B_I_eq = B(self.v_idx{i}, self.I{i});
                    B_I_ineq = B(~self.v_idx{i}, self.I{i});

                    AB_I = [A_I; B_I_eq];
                    [AB_I, self.D_idx{i}] = linearly_independent_rows(AB_I);
                    [m_I, n_I] = size(AB_I);
                    if rank(AB_I) ~= m_I
                        error('Something is wrong. Matrix should have had linearly independent rows');
                    end

                    if n_I == m_I
                        self.A_subcase{i} = 1;
                        self.D{i} = inv(AB_I);
                    elseif n_I == m_I + 1 && isempty(B_I_ineq)
                        % TODO: implement for B
                        self.A_subcase{i} = 2;
                        self.D{i} = AB_I' / (AB_I*AB_I');
                        self.a{i} = null(AB_I);
                    else
                        self.A_subcase{i} = 3;
                        self.pars_subset{i} = Pars(AB_I, B_I_ineq, remove_same_columns, self.tol);
                    end
                    AB = [A; B];
                    self.d_vec{i} = sum(AB(:, self.J{i}), 2) - sum(AB(:, self.K{i}), 2);
                end
            end
        end

        function check_A(self)
            if self.m_A == 0
                error('Matrix A must not be empty');
            end
            if rank(self.A_original) ~= self.m_A
                error('Matrix does not have linearly independent rows');
            end
        end

        function [A, B] = modify_input_matrix(self, A, B, remove_same_columns)
            self.A_original = A;
            self.B_original = B;
            [self.m_A, self.n] = size(A);
            [self.m_B, ~] = size(B);

            % Find zeros columns
            self.zero_columns = vecnorm([A; B], 2, 1) <= self.tol;
            A = A(:, ~self.zero_columns);
            B = B(:, ~self.zero_columns);

            if remove_same_columns
                % Find columns which are multiples of each other
                [self.multiples, self.multiples_counts] = find_column_multiples([A; B], self.tol);
                % Multiply the columns which are multiplied
                for i = 1:size(self.multiples_counts,1)
                    k = self.multiples_counts(i, 1);
                    A(:,k) = A(:,k) * self.multiples_counts(i, 2);
                    B(:,k) = B(:,k) * self.multiples_counts(i, 2);
                end
                % Remove columns which are multiples
                A(:, self.multiples(:,2)) = [];
                B(:, self.multiples(:,2)) = [];
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
