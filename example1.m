clear all;
close all;

addpath("src")

%% Input: Clarke's transform A

A = 2/5*[1, cos(2/5*pi), cos(4/5*pi), cos(-4/5*pi), cos(-2/5*pi);...
         0, sin(2/5*pi), sin(4/5*pi), sin(-4/5*pi), sin(-2/5*pi)];    

%% Input: Required voltage y

ts = 0:100e-6:0.04;

[n_y, n_x] = size(A);
n_t = length(ts);

ys = zeros(n_y, n_t);
for k = 1:n_t
    wt = 2*pi*50*ts(k);
    ys(:,k) = 1*[cos(wt); sin(wt)];
end

%% Output: Precompute offline part

pars = Pars(A);
solver = Solver(pars);

%% Output: Get a solution for every time

xs = zeros(n_x, n_t);
xs_l2 = zeros(n_x, n_t);
A_inv = A' / (A*A');
for k = 1:n_t
    xs(:,k) = solver.min_effort(ys(:,k));
    xs_l2(:,k) = A_inv*ys(:,k);
end

%% Plot the results

ylims = [-max(xs_l2(:)); max(xs_l2(:))];

fig = figure();
plot(ts, ys');
xlabel('Time [s]');
ylabel('Required voltage');
saveas(fig, 'figures/res1.png')

fig = figure();
plot(ts, xs');
xlabel('Time [s]');
ylabel('Input voltage');
title('Our approach');
ylim(ylims);
saveas(fig, 'figures/res2.png')

fig = figure();
plot(ts, xs_l2');
xlabel('Time [s]');
ylabel('Input voltage');
title('Standard approach');
ylim(ylims);
saveas(fig, 'figures/res3.png')
