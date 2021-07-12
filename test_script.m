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
A = get_a(A_type);




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
