function U = get_u(A)
    % get_u Computes the set of extremal points of the dual feasible set
    %    ||A'*u|| <= 1
    %
    % Inputs:
    % A (matrix): matrix A specifying the problem.
    %
    % Outputs:
    % U (matrix): extremal points of the dual feasible set.
    
    tol = 1e-10;    
    [m, n] = size(A);

    % Constraints in the enhanced space B*[u1^+ u1^- u2tilde z^+ z^-] = b
    C = [A' -A' -eye(n) eye(n); ...
        zeros(1,m) zeros(1,m) ones(1,n) ones(1,n)];
    c = [zeros(n,1); 1];

    % Generate all indices to select square submatrices of B
    iis = nchoosek(1:size(C,2),size(C,1));

    % Run over all submatrices and get the extremal points in the enhanced space
    U_ext = [];
    for i = 1:size(iis,1)
        U_ext = [U_ext, Get_Corner(C, c, iis(i,:))];
    end

    % Reduce the extremal points to the original space
    U = [U_ext(1:m,:) - U_ext(m+1:2*m,:)];

    % Remove duplicities
    U = unique(U', 'rows')';

    % Remove points which are not extremal in the original space
    ii = convhulln(U');
    ii = unique(ii(:));
    U = U(:,ii)';

    % This is used only for numerical errors
    U = convhulln_error(U);

    % Feasibility check for ||A'*u|| = 1
    if max(abs(sum(abs(A'*U')) - 1)) >= tol
        error('Solutions do not have norm 1');
    end
end


function u = Get_Corner(C, c, ii)
    % Get_Corner computes the extremal set of the set Bv = b, v >= 0.

    tol = 1e-10;
    u = [];
    B_red = C(:,ii);
    if rank(B_red) == size(B_red,1)
        u_red = B_red \ c;
        if min(u_red) >= -tol
            u_red(u_red <= tol) = 0;
            u = zeros(size(C,2),1);
            u(ii) = u_red;
        end
    end
end

function U = convhulln_error(U)
    % convhulln_error returns the convex hull while allowing for epsilon errors.

    tol = 1e-10;
    remove = false(size(U,1),1);
    for i = 1:size(U,1)
        u = U(i,:);
        d = U - u;
        d = d ./ sum(abs(d), 2);
        [d_sort, ii_sort] = sortrows(d);

        for j = 1:size(d_sort,1)-1
            if norm(d_sort(j,:) - d_sort(j+1,:)) <= tol
                j1 = ii_sort(j);
                j2 = ii_sort(j+1);
                d1 = U(j1,:) - u;
                d2 = U(j2,:) - u;
                d1_norm = sum(abs(d1));
                d2_norm = sum(abs(d2));

                if norm(d1/d1_norm - d2/d2_norm) > tol
                    error("Something is wrong");
                end

                if d1_norm > d2_norm
                    remove(j2) = true;
                else
                    remove(j1) = true;
                end
            end
        end
    end

    U = U(~remove,:);
end
