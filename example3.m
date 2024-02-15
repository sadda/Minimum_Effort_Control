% clear all;
% close all;

addpath("src")
addpath("user-provided")

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

t = 0;
dt = 20e-6;
Tsim = 0.06;
N = round(Tsim/dt) + 1;

% Compute or load U
file_name = 'data/pars_3.mat';
if ~isfile(file_name)
    pars = get_u(A);
    if ~isfolder('data')
        mkdir('data');
    end
    save(file_name, "pars");
else
    load(file_name, "pars");
end

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

    % [x, ~, pars] = min_effort(pars, y, @find_x_3);
    [x, ~, pars] = min_effort(pars, y, []);

    o_t(k,1) = t;
    o_y(k,:) = y';
    o_x(k,:) = x';
    t = t + dt;
end

y_lim = 0.8;

figure;
subplot(2,2,1);
plot(o_t, o_y, 'LineWidth', 3); grid on;

subplot(2,2,2);
plot(o_t, [o_x(:,1),o_x(:,2),o_x(:,3)], 'LineWidth', 3); grid on;
ylim([-y_lim,y_lim]);

subplot(2,2,3);
plot(o_t, [o_x(:,4),o_x(:,5),o_x(:,6)], 'LineWidth', 3); grid on;
ylim([-y_lim,y_lim]);

subplot(2,2,4);
plot(o_t, [o_x(:,7),o_x(:,8),o_x(:,9)], 'LineWidth', 3); grid on;
ylim([-y_lim,y_lim]);

