clear all;
close all;

for A_type = 1:7
    
    fprintf('Running test: %d\n', A_type);
    
    %% Input: Clarke's transform A
    
    A = get_a(A_type);
    [m, n] = size(A);
    
    %% Input: Required voltage y
    
    ts = 0:100e-6:0.06;
    N = length(ts);
    
    omega = 2*pi*50;
    Um = 230*sqrt(2);
    
    ys = zeros(m,N);
    for k = 1:N
        wt = omega*ts(k);
        
        ys(1,k) = Um*cos(wt);
        ys(2,k) = Um*sin(wt);
        if A_type == 6
            ys(3,k) = -0.4*Um*cos(wt);
        end
    end
    
    %% Output: Precompute set U
    
    U = get_u(A);
    
    %% Output: Get a solution for every time
    
    xs1 = zeros(n,N);
    xs2 = zeros(n,N);
    for k = 1:N
        xs1(:,k) = min_effort(A, ys(:,k), U);
        xs2(:,k) = min_effort_linprog(A, ys(:,k));
    end
    
    %% Check for equality
    
    if norm(xs1 - xs2) <= 1e-10
        fprintf('Test passed. Difference = %1.3e\n', norm(xs1 - xs2));
    else
        error('Test failed for A_type = %d. Difference = %1.3e.', A_type, norm(xs1 - xs2))
    end
end
