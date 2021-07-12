function corners = get_u(A)
    % get_u computes the set of extremal points of the dual feasible set ||A'*u|| <= 1.
    %
    % Inputs:
    % A (matrix): matrix A specifying the problem.
    %
    % Outputs:
    % U (matrix): extremal points of the dual feasible set.
    
    [m, n] = size(A);
    
    % Constraints in the enhanced space B*[u^+ u^- z^+ z^- w] = b
    B = [A' -A' -eye(n) eye(n) zeros(n,1); ...
        zeros(1,m) zeros(1,m) ones(1,n) ones(1,n) 1];
    b = [zeros(n,1); 1];
    
    % Generate all indices to select square submatrices of B
    iis = nchoosek(1:size(B,2),size(B,1));
    
    % Run over all submatrices and get the extremal points in the enhanced space
    corners_ext = [];
    for i = 1:size(iis,1)
        corners_ext = [corners_ext, Get_Corner(B, b, iis(i,:))];
    end
    
    % Reduce the extremal points to the original space
    corners = corners_ext(1:m,:) - corners_ext(m+1:2*m,:);
    
    % Remove duplicities
    corners = unique(corners', 'rows')';
    
    % Remove points which are not extremal in the original space
    ii = convhulln(corners');
    ii = unique(ii(:));
    corners = corners(:,ii)';
    
    % This is used only for numerical errors
    corners = convhulln_error(corners);
end


function u = Get_Corner(B, b, ii)
    % Get_Corner computes the extremal set of the set Bv = b, v >= 0.
    
    tol = 1e-10;
    u = [];
    B_red = B(:,ii);
    if rank(B_red) == size(B_red,1)
        u_red = B_red \ b;
        if min(u_red) >= -tol
            u_red(u_red <= tol) = 0;
            u = zeros(size(B,2),1);
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


