clear all;
close all;

addpath("src")

f = 50;
w = 2*pi*f;
fi = 2/3*pi;

A1 = 2/3*[1, cos(fi), cos(-fi);...
    0, sin(fi), sin(-fi);...
    1/2, 1/2, 1/2];

% 3x trifazovy stridac
% menic 1: u11, u12, u13
% menic 2: u21, u22, u23
% menic 3: u31, u32, u33

%x = [u11, u12, u13, u21, u22, u23, u31, u32, u33]'
A2 = [ 1  0 -1  0 -1  1  0  0  0;... % ua = u11 - u13 + u23 - u22
    0  0  0  1  0 -1  0 -1  1;... % ub = u21 - u23 + u33 - u32
    0 -1  1  0  0  0  1  0 -1;];  % uc = u31 - u33 + u13 - u12

A = A1*A2;
B = [ones(1, 9); -ones(1, 9)];

% A = [
%     0 1 2 3 4 5 6 7 8;
%     2 3 4 5 6 7 8 9 9;
% ];
% A = randn(4,8);
% A = randn(2,2);
% A = [A A A];

% Compute or load U

% A(:,[5,7,8]) = [];
pars = Pars(A, B, true);
solver = Solver(pars);


t = 0;
dt = 20e-6;
Tsim = 0.06;
N = round(Tsim/dt) + 1;

for k = 1:N
    wt = w*t;
    % Three phase voltage, ground fault (zero-sequence) in t = ....
    if t<0.025
        ua = 1*cos(wt);
        ub = 1*cos(wt - 2/3*pi);
        uc = 1*cos(wt + 2/3*pi);
    else
        ua = 1*sqrt(3)*cos(wt);
        ub = 1*sqrt(3)*cos(wt-pi/3);
        uc = 0;
    end
    y = A1*[ua;ub;uc];
    %y = [y; 10000*ones(size(B,1),1)];
    y = [y; 0.3*ones(size(B,1),1)];

%     y = [y;0];
% y = y(1:size(A,1));


    x = solver.min_effort(y);
%     x = solver.min_effort_user_provided(y);

    o_t(k,1) = t;
    o_y(k,:) = y';
    o_x(k,:) = x';
    t = t + dt;
end

% y_lim = 0.8;

fig = figure;
plot(o_t, [o_x(:,1),o_x(:,2),o_x(:,3)], 'LineWidth', 3); grid on;
% ylim([-y_lim,y_lim]);
xlim([0,Tsim]);
set(gcf, 'Position', get(0, 'Screensize'));
% exportgraphics(fig, 'Res3.png', 'Resolution', 600)
% 
% pars.plot_s_min_s_max(dt, 'LineWidth', 3)
% xlim([0,Tsim]);
% set(gcf, 'Position', get(0, 'Screensize'));
% 
% pars.analyze_solution(true);
% 
