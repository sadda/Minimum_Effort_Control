% Example of usage. See the GitHub readme for more comments.
%
% https://github.com/sadda/Minimum_Effort_Control

clear all;
close all;

%% Settings

plot_figure = true;
kappas = linspace(0.4, 0.01, 40);

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

%% Input: Required voltage y

ts = 0:100e-6:0.06;
omega = 2*pi*50;
Um = 230*sqrt(2);

[n_y, n_x] = size(A);
n_t = length(ts);
n_k = length(kappas);

ys = zeros(n_y, n_t, n_k);
for i = 1:length(kappas)
    for k = 1:n_t
        wt = omega*ts(k);
        
        ys(1,k,i) = Um*cos(wt);
        ys(2,k,i) = Um*sin(wt);
        ys(3,k,i) = abs(kappas(i)*Um*cos(3*wt));
        ys(4,k,i) = abs(kappas(i)*Um*sin(3*wt));
        ys(5,k,i) = ys(3,k,i);
        ys(6,k,i) = ys(4,k,i);
    end
end

%% Output: Precompute set U

U = get_u(A1, A2);

%% Output: Get a solution for every time

xs = zeros(n_x, n_t, n_k);
for i = 1:length(kappas)
    for k = 1:n_t
        xs(:,k,i) = min_effort(A1, A2, ys(:,k,i), U);
    end
end

%% Plot the results

if plot_figure
    for i = 1:length(kappas)
        y_lim = [-400 400];
        h = plot_figures(ts, xs(:,:,i), A1, A2, Um, y_lim, kappas(i));
        
        filename = 'Result.gif';
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        if i == 1
            imwrite(imind,cm,filename,'gif','Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
        end
    end
end


