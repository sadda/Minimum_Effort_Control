% Example of usage. See the GitHub readme for more comments.
%
% https://github.com/sadda/Minimum_Effort_Control

clear all;
close all;

%% Settings

plot_figure = true;
% kappa_all = linspace(0.9, 0.01, 90);
kappa_all = linspace(0.5, 0.02, 30);
omega_coef = 3;

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
    kappa
    ys = zeros(m,N);
    for k = 1:N
        wt = omega*ts(k);
        
        ys(1,k) = Um*cos(wt);
        ys(2,k) = Um*sin(wt);
        ys(3,k) = abs(kappa*Um*cos(omega_coef*wt));
        ys(4,k) = abs(kappa*Um*sin(omega_coef*wt));
        ys(5,k) = ys(3,k);
        ys(6,k) = ys(4,k);
    end
    
    %% Output: Get a solution for every time
    
    xs = zeros(n,N);
    for k = 1:N
        xs(:,k) = min_effort(A1, A2, ys(:,k), U);
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
        plot(Um*kappa*cos(disc), Um*kappa*sin(disc), 'k--');
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



