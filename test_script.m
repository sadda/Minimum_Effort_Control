clear all; close all;
% Test script to test functions provided by Lukas Adam
f = 50;
w = 2*pi*f;

dt = 100e-6;
t = 0;
T = 0.06;

N = round(T/dt);

Um = 1; %230*sqrt(2);

A_type = 6;

switch A_type
    case 1
        % 3-faz
        A = 2/3*[1, cos(2/3*pi), cos(-2/3*pi);...
            0, sin(2/3*pi), sin(-2/3*pi);];
    case 2
        % 5-faz 3 DOF
        A = 2/5*[1, cos(2*pi/5), cos(4*pi/5), cos(-4*pi/5), cos(-2*pi/5); ...
            0, sin(2*pi/5), sin(4*pi/5), sin(-4*pi/5), sin(-2*pi/5);];
    case 3
        % 5-faz 1 DOF
        A = 2/5*[1, cos(2*pi/5), cos(4*pi/5), cos(-4*pi/5), cos(-2*pi/5); ...
            0, sin(2*pi/5), sin(4*pi/5), sin(-4*pi/5), sin(-2*pi/5); ...
            1, cos(4*pi/5), cos(-2*pi/5), cos(2*pi/5), cos(-4*pi/5); ...
            0, sin(4*pi/5), sin(-2*pi/5), sin(2*pi/5), sin(-4*pi/5)];
    case 4
        % 7-faz 5 DOF
        alpha = 2*pi/7;
        A = 2/7*[1, cos(alpha), cos(2*alpha), cos(3*alpha), cos(4*alpha), cos(5*alpha), cos(6*alpha); ...
            0, sin(alpha), sin(2*alpha), sin(3*alpha), sin(4*alpha), sin(5*alpha), sin(6*alpha);];
    case 5
        % 7-faz 1 DOF
        alpha = 2*pi/7;
        A = 2/7*[1, cos(alpha), cos(2*alpha), cos(3*alpha), cos(4*alpha), cos(5*alpha), cos(6*alpha);...
            0, sin(alpha), sin(2*alpha), sin(3*alpha), sin(4*alpha), sin(5*alpha), sin(6*alpha);...
            1, cos(2*alpha), cos(4*alpha), cos(6*alpha), cos(8*alpha), cos(10*alpha), cos(12*alpha);...
            0, sin(2*alpha), sin(4*alpha), sin(6*alpha), sin(8*alpha), sin(10*alpha), sin(12*alpha);...
            1, cos(3*alpha), cos(6*alpha), cos(9*alpha), cos(12*alpha), cos(15*alpha), cos(18*alpha);...
            0, sin(3*alpha), sin(6*alpha), sin(9*alpha), sin(12*alpha), sin(15*alpha), sin(18*alpha);];
    case 6
        % 4-faz 1 DOF
        A1 = 2/3*[1, cos(2/3*pi), cos(-2/3*pi);...
            0, sin(2/3*pi), sin(-2/3*pi);...
            1/2, 1/2, 1/2;];
        
        A2 = [1, 0, 0, -1;...
            0, 1, 0, -1;...
            0, 0, 1, -1;];
        A = A1*A2;
        
end



[m, n] = size(A);
y = zeros(m,1);

U = get_u(A);


for k = 1:N
    wt = w*t;
    
    Uaf = Um*cos(wt);
    Ubt = Um*sin(wt);
    U0 = -0.4*cos(wt);
    
    y(1) = Uaf;
    y(2) = Ubt;
    if A_type == 6
        y(3) = U0;
    end
    
    % Solution found by linear programming
    x_lp = solve_linprog(A, y);
    
    % Solution found by Lukas Adam's method
    x_ad = solve_ours(A, y, U);
    
    o_t(k,1) = t;
    o_y(k,:) = y';
    o_xlp(k,:) = x_lp';
    o_xad(k,:) = x_ad';
    t = t + dt;
end


figure;
subplot(3,1,1);
plot(o_t, o_y); grid on;
title('y(1),y(2),y(3)');

subplot(3,1,2);
plot(o_t, o_xlp); grid on;
title('x(1),x(2),x(3),x(4) by linear programming');

subplot(3,1,3);
plot(o_t, o_xad); grid on;
title('x(1),x(2),x(3),x(4) by ADMM');

% No more
