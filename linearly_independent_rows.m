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
