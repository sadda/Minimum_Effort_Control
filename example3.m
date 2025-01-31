clear all;
close all;

addpath("src")

f = 50;
w = 2*pi*f;
fi = 2/3*pi;

A1 = 2/3*[1, cos(fi), cos(-fi);...
    0, sin(fi), sin(-fi);...
    1/2, 1/2, 1/2];
A2 = [ 1  0 -1  0 -1  1  0  0  0;... % ua = u11 - u13 + u23 - u22
    0  0  0  1  0 -1  0 -1  1;... % ub = u21 - u23 + u33 - u32
    0 -1  1  0  0  0  1  0 -1;];  % uc = u31 - u33 + u13 - u12

A = A1*A2;

t = 0;
dt = 20e-6;
Tsim = 0.06;
N = round(Tsim/dt) + 1;

% Precompute the offline part
pars = Pars(A);
solver = Solver(pars);

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
    x = solver.min_effort(y);

    o_t(k) = t;
    o_y(k,:) = y';
    o_x(k,:) = x';
    t = t + dt;
end

y_lim = 0.8;

fig = figure;
plot(o_t, [o_x(:,1),o_x(:,2),o_x(:,3)], 'LineWidth', 3); grid on;
ylim([-y_lim,y_lim]);
xlim([0,Tsim]);
set(gcf, 'Position', get(0, 'Screensize'));
