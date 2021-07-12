clear all;
close all;

%% Input: Clarke's transform A

A_type = 6;
A = get_a(A_type);
[m, n] = size(A);

%% Input: Required voltage y

ts = 0:100e-6:0.06;
N = length(ts);

omega = 2*pi*50;
Um = 230*sqrt(2);

ys = zeros(m,N);
for k = 1:N
    wt = omega*ts(k);
    
    ys(1,k) = Um*cos(wt);
    ys(2,k) = Um*sin(wt);
    if A_type == 6
        ys(3,k) = -0.4*Um*cos(wt);
    end
end

%% Output: Precompute set U

U = get_u(A);

%% Output: Get a solution for every time

xs1 = zeros(n,N);
for k = 1:N
    xs1(:,k) = solve_ours(A, ys(:,k), U);
end

%% Output: Get a solution for every time

xs2 = zeros(n,N);
for k = 1:N
    xs2(:,k) = A'*((A*A') \ ys(:,k));
end

%% Output: Plot the results

ylims = [min([xs1(:); xs2(:)]); max([xs1(:); xs2(:)])];

figure();
plot(ts, ys');
grid on;
xlabel('Time');
ylabel('Required voltage');
print(gcf, 'res1.png', '-dpng', '-r600');

figure();
plot(ts, xs1');
grid on;
xlabel('Time');
ylabel('Input voltage');
title('Our approach');
ylim(ylims);
print(gcf, 'res2.png', '-dpng', '-r600');

figure();
plot(ts, xs2');
grid on;
xlabel('Time');
ylabel('Input voltage');
title('Standard approach');
ylim(ylims);
print(gcf, 'res3.png', '-dpng', '-r600');
