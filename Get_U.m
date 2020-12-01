function corners = Get_U(A)
    
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
end


function u = Get_Corner(B, b, ii)
    
    u = [];
    B_red = B(:,ii);
    if rank(B_red) == size(B_red,1)
        u_red = B_red \ b;
        if min(u_red) >= -1e-10
            u_red(u_red <= 1e-10) = 0;
            u = zeros(size(B,2),1);
            u(ii) = u_red;
        end
    end
end

