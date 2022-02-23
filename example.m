% Example of usage. See the GitHub readme for more comments.
%
% https://github.com/sadda/Minimum_Effort_Control

clear all;
close all;

%% Input: Clarke's transform A

A_type = 3;
A = get_a(A_type);

i1 = [1 2];
i2 = [3 4];
A1 = A(i1,:);
A2 = [A(i2,:); -A(i2,:)];

A = [A1; A2];
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
%     if A_type == 6
%         ys(3,k) = -0.4*Um*cos(wt);
%     end
    ys(3,k) = 0.1+abs(0.1*Um*cos(wt));
    ys(4,k) = 0.1+abs(0.1*Um*sin(wt));
    ys(5,k) = ys(3,k);
    ys(6,k) = ys(4,k);
end

%% Output: Precompute set U

U = get_u(A1, A2);

%% Output: Get a solution for every time

xs1 = zeros(n,N);
for k = 1:N
    xs1(:,k) = min_effort(A1, A2, ys(:,k), U);
end

%% Output: Get the standard solution

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
% print(gcf, 'res1.png', '-dpng', '-r600');

figure();
plot(ts, xs1');
grid on;
xlabel('Time');
ylabel('Input voltage');
title('Our approach');
ylim(ylims);
% print(gcf, 'res2.png', '-dpng', '-r600');

figure();
plot(ts, xs2');
grid on;
xlabel('Time');
ylabel('Input voltage');
title('Standard approach');
ylim(ylims);
% print(gcf, 'res3.png', '-dpng', '-r600');


figure();
plot(ts, (A1*xs1)');


figure();
plot(ts, (A2(1:2,:)*xs1)');

