clear all;

A_type = 1;

switch A_type
    case 1
        A = 2/3*[1, cos(2*pi/3) cos(-2*pi/3); ...
            0, sin(2*pi/3), sin(-2*pi/3)];
    case 2
        A = 2/5*[1, cos(2*pi/5), cos(4*pi/5), cos(-4*pi/5), cos(-2*pi/5); ...
            0, sin(2*pi/5), sin(4*pi/5), sin(-4*pi/5), sin(-2*pi/5); ...
            1, cos(4*pi/5), cos(-2*pi/5), cos(2*pi/5), cos(-4*pi/5); ...
            0, sin(4*pi/5), sin(-2*pi/5), sin(2*pi/5), sin(-4*pi/5)];
    case 3
        A = randn(4,5);
    case 4
        A = randn(7,8);
    case 5
        A = 2/5*[1, cos(2*pi/5), cos(4*pi/5), cos(-4*pi/5), cos(-2*pi/5); ...
            0, sin(2*pi/5), sin(4*pi/5), sin(-4*pi/5), sin(-2*pi/5)];
    case 6
        A = randn(3,8);
end

mn = size(A);
m  = mn(1);
n  = mn(2);


% y = rand(m,1);
y = [-325.2691; 0];

[x1, val1] = Solve_Linprog(A, y, 1);
[x2, val2] = Solve_Linprog(A, y, 0);
[x3, val3] = Solve_Ours1(A, y);
[x4, val4] = Solve_Ours2(A, y);
[x5, val5] = Solve_Ours3(A, y);

[x1, x2, x3, x4, x5]





