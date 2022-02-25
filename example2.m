% Example of usage. See the GitHub readme for more comments.
%
% https://github.com/sadda/Minimum_Effort_Control

% clear all;
close all;

%% Settings

box = false;
plot_figure = true;
kappa_all = linspace(0.9, 0.01, 90);
omega_coef = 3;
angle = 0;

%% Input: Clarke's transform A

phi5 = 2/5*pi;
A_orig = 2/5*[1, cos(phi5), cos(2*phi5), cos(-2*phi5), cos(-phi5); ...
    0, sin(phi5), sin(2*phi5), sin(-2*phi5), sin(-phi5); ...
    1, cos(2*phi5), cos(-phi5), cos(phi5), cos(-2*phi5); ...
    0, sin(2*phi5), sin(-phi5), sin(phi5), sin(-2*phi5)];

i1 = [1 2];
i2 = [3 4];
A1 = A_orig(i1,:);
A2 = [A_orig(i2,:); -A_orig(i2,:)];
A = [A1; A2];

[m, n] = size(A);

%% Output: Precompute set U

U = get_u(A1, A2);

%% Input: Required voltage y

ts = 0:100e-6:0.06;
N = length(ts);

omega = 2*pi*50;
Um = 230*sqrt(2);

for kappa = kappa_all
    ys = zeros(m,N);        
    for k = 1:N
        wt = omega*ts(k);
        
        ys(1,k) = Um*cos(wt);
        ys(2,k) = Um*sin(wt);
        if box
            ys(3,k) = kappa*Um;
            ys(4,k) = kappa*Um;
        else
            ys(3,k) = abs(kappa*Um*cos(omega_coef*wt+angle));
            ys(4,k) = abs(kappa*Um*sin(omega_coef*wt+angle));
        end
        ys(5,k) = ys(3,k);
        ys(6,k) = ys(4,k);
    end
    
    %% Output: Get a solution for every time
    
    xs = zeros(n,N);
    for k = 1:N
        %     xs(:,k) = min_effort(A1, A2, ys(:,k), U);
        
        tol = 1e-10;
        y = ys(:,k);
        
        % Compute the optimal dual solution
        [val, i_max] = max(U*y);
        u_opt = U(i_max, :)';
        
        % Assign the index sets for complementarity
        I0 = abs(A'*u_opt) <= tol;
        I1 = A'*u_opt > tol;
        I2 = A'*u_opt < -tol;
        J = [true(size(A1,1), 1); u_opt(size(A1,1)+1:end) < -tol];
        
        % Use the complementarity conditions to compute the primal solution
        x = zeros(n,1);
        x(I1) = val;
        x(I2) = -val;
        
        % Hard-coded bad conditioning
        if all(I0 == logical([1;0;1;1;0])) && (all(J == logical([1;1;0;0;0;1])) || all(J == logical([1;1;0;1;0;0])))
            t_p = A([1 2],I0) \ (y([1 2]) - A([1 2],I1|I2)*x(I1|I2));
            t_d = null(A(J,I0));
            if any(t_d <= 0)
                if all(t_d <= 0)
                    t_d = -t_d;
                else
                    error('Not yet implemented.');
                end
            end
            if A(3,I0)*t_d <= 0
                error('Not yet implemented.');
            end
            c_max1 = min((val-t_p) ./ t_d);
            c_min1 = max((-val-t_p) ./ t_d);
            c_max2 = (y(3) - A(3,I1|I2)*x(I1|I2) - A(3,I0)*t_p) / (A(3,I0)*t_d);
            c_min2 = (-y(3) - A(3,I1|I2)*x(I1|I2) - A(3,I0)*t_p) / (A(3,I0)*t_d);
            c_max = min(c_max1, c_max2);
            c_min = max(c_min1, c_min2);
            if c_max <= c_min - tol
                error('Something wrong');
            end
            
            xs_a = x;
            xs_a(I0) = t_p + c_min*t_d;
            xs_b = x;
            xs_b(I0) = t_p + c_max*t_d;
            
            norm_a = norm(A([3 4],:)*xs_a);
            norm_b = norm(A([3 4],:)*xs_b);
                        
            if norm(norm_a) > norm(norm_b) + tol
                c = c_min;
            elseif norm(norm_b) > norm(norm_a) + tol
                c = c_max;
            elseif k > 1 && norm(xs_a - xs(:,k-1)) <= norm(xs_b - xs(:,k-1)) - tol
                c = c_min;
            elseif k > 1 && norm(xs_b - xs(:,k-1)) <= norm(xs_a - xs(:,k-1)) - tol
                c = c_max;
            else
                c = 0.5*(c_min+c_max);
            end
            
            x(I0) = t_p + c*t_d;
        else
            x(I0) = A(J,I0) \ (y(J) - A(J,I1|I2)*x(I1|I2));
        end
        xs(:,k) = x;
    end
    
    if plot_figure
        disc = linspace(0, 2*pi, 100);
        y_lim =[-400 400];
        
        h = figure('visible', 'off', 'position', 4*[0, 0, 400, 100]);
        axis tight manual % this ensures that getframe() returns a consistent size
        
        subplot(1,4,1);
        plot(ts, xs');
        ylim(y_lim);
        ylabel('kappa = ' + string(kappa))
        title('xs');
        
        subplot(1,4,2);
        plot(ts, (A1*xs)');
        ylim(y_lim);
        title('A*xs');
        
        subplot(1,4,3);
        hold on;
        plot(ts, (A2([1 2],:)*xs)');
        plot(ts, Um*kappa*ones(1,length(ts)), 'k--');
        plot(ts, -Um*kappa*ones(1,length(ts)), 'k--');
        ylim(y_lim);
        title('B*xs');
        
        subplot(1,4,4);
        hold on;
        plot((A1(1,:)*xs)', (A1(2,:)*xs)');
        plot((A2(1,:)*xs)', (A2(2,:)*xs)');
        if box
            plot([-Um*kappa; Um*kappa; Um*kappa; -Um*kappa; -Um*kappa], [-Um*kappa; -Um*kappa; Um*kappa; Um*kappa; -Um*kappa], 'k--')
        else
            plot(Um*kappa*cos(disc), Um*kappa*sin(disc), 'k--');
        end
        xlim(y_lim);
        ylim(y_lim);
        title('B*xs');
        
        drawnow
        
        filename = 'Result.gif';
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        if kappa == kappa_all(1)
            imwrite(imind,cm,filename,'gif','Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
        end
    end
end



