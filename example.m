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
% U = get_u([A1; A2(1:2,:)]);

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
        %         ys(3,k) = abs(kappa*Um*cos(omega_coef*wt));
        %         ys(4,k) = abs(kappa*Um*sin(omega_coef*wt));
        ys(3,k) = kappa*Um;
        ys(4,k) = kappa*Um;
        ys(5,k) = ys(3,k);
        ys(6,k) = ys(4,k);
    end
    
    %% Output: Get a solution for every time
    
    xs = zeros(n,N);
    for k = 1:N
        xs(:,k) = min_effort(A1, A2, ys(:,k), U);
        %         xs(:,k) = min_effort([A1; A2(1:2,:)], [], ys(1:4,k), U);
    end
    
    if plot_figure
        y_lim = [-400 400];
        h = plot_figures(ts, xs, A1, A2, Um, y_lim, kappa);
        
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



