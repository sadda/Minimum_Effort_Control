clc;
clear all;
close all;

addpath("src")

%% Time
ts = 0:100e-6:0.04;
n_t = length(ts);
wt = 2*pi*50*ts;
remove_same_columns = true;
y_original = true;

%% Input: Clarke's transform A and required voltage y

for j = 1:18
    j
    switch j
        case 1
            % Three-phase system, DOF: V0
            fi = 2/3*pi;
            A = 2/3*[1, cos(fi), cos(-fi);...
                0, sin(fi), sin(-fi)];
            ys = 1*[cos(wt); sin(wt)];
        case 2
            % Three-phase four-leg system, DOF: V0
            fi = 2/3*pi;
            A1 = 2/3*[1, cos(fi), cos(-fi);...
                0, sin(fi), sin(-fi)];
            A2 = [1, 0, 0, -1;...
                0, 1, 0, -1;...
                0, 0, 1, -1];
            A = A1*A2;
            ys = 1*[cos(wt); sin(wt)];
        case 3
            % Three-phase four-leg system, defined V0
            fi = 2/3*pi;
            A1 = 2/3*[1, cos(fi), cos(-fi);...
                0, sin(fi), sin(-fi);...
                1/2, 1/2, 1/2];
            A2 = [1, 0, 0, -1;...
                0, 1, 0, -1;...
                0, 0, 1, -1];
            A = A1*A2;
            ys = 1*[cos(wt); sin(wt); -cos(wt)];
        case 4
            % Five-phase system, DOF: Vx1, Vy1, V0
            fi = 2/5*pi;
            A = 2/5*[1, cos(fi), cos(2*fi), cos(-2*fi), cos(-fi);...
                0, sin(fi), sin(2*fi), sin(-2*fi), sin(-fi)];
            ys = 1*[cos(wt); sin(wt)];
        case 5
            % Five-phase system, DOF: V0
            fi = 2/5*pi;
            A = 2/5*[1, cos(fi), cos(2*fi), cos(-2*fi), cos(-fi);...
                0, sin(fi), sin(2*fi), sin(-2*fi), sin(-fi);...
                1, cos(-2*fi), cos(fi), cos(-fi), cos(2*fi);...
                0, sin(-2*fi), sin(fi), sin(-fi), sin(2*fi)];
            if y_original
                ys = 1*[cos(wt); sin(wt); 0.15*cos(3*wt); 0.15*sin(3*wt)];
            else
                ys = 1*[cos(wt); sin(wt); zeros(1,n_t); zeros(1,n_t)];
            end
        case 6
            % Five-phase system, fault-tolerant mode, DOF: Vx1, Vy1, V0
            fi = 2/5*pi;
            A = 2/5*[0, cos(fi), cos(2*fi), cos(-2*fi), cos(-fi);...
                0, sin(fi), sin(2*fi), sin(-2*fi), sin(-fi)];
            ys = 1*[cos(wt); sin(wt)];
        case 7
            % Five-phase system, multilevel converter, fault-tolerant mode, DOF: Vx1, Vy1, V0
            fi = 2/5*pi;
            A = 2/5*[0.5, cos(fi), cos(2*fi), cos(-2*fi), cos(-fi);...
                0, sin(fi), sin(2*fi), sin(-2*fi), sin(-fi)];
            ys = 1*[cos(wt); sin(wt)];
        case 8
            % Seven-phase system, DOF: Vx1, Vy1, Vx2, Vy2, V0
            fi = 2/7*pi;
            A = 2/7*[1, cos(fi), cos(2*fi), cos(3*fi), cos(-3*fi), cos(-2*fi), cos(-fi);...
                0, sin(fi), sin(2*fi), sin(3*fi), sin(-3*fi), sin(-2*fi), sin(-fi);];

            ys = 1*[cos(wt); sin(wt)];
        case 9
            % Seven-phase system, DOF: Vx2, Vy2, V0
            fi = 2/7*pi;
            A = 2/7*[1, cos(fi), cos(2*fi), cos(3*fi), cos(-3*fi), cos(-2*fi), cos(-fi);...
                0, sin(fi), sin(2*fi), sin(3*fi), sin(-3*fi), sin(-2*fi), sin(-fi);...
                1, cos(3*fi), cos(-fi), cos(2*fi), cos(-2*fi), cos(fi), cos(-3*fi);...
                0, sin(3*fi), sin(-fi), sin(2*fi), sin(-2*fi), sin(fi), sin(-3*fi)];
            if y_original
                ys = 1*[cos(wt); sin(wt); 0.15*cos(3*wt); 0.15*sin(3*wt)];
            else
                ys = 1*[cos(wt); sin(wt); zeros(1,n_t); zeros(1,n_t)];
            end
        case 10
            % Seven-phase system, DOF: Vx1, Vy1, V0
            fi = 2/7*pi;
            A = 2/7*[1, cos(fi), cos(2*fi), cos(3*fi), cos(-3*fi), cos(-2*fi), cos(-fi);...
                0, sin(fi), sin(2*fi), sin(3*fi), sin(-3*fi), sin(-2*fi), sin(-fi);...
                1, cos(-2*fi), cos(3*fi), cos(fi), cos(-fi), cos(-3*fi), cos(2*fi);...
                0, sin(-2*fi), sin(3*fi), sin(fi), sin(-fi), sin(-3*fi), sin(2*fi)];
            if y_original
                ys = 1*[cos(wt); sin(wt); 0.08*cos(5*wt); 0.08*sin(5*wt)];
            else
                ys = 1*[cos(wt); sin(wt); zeros(1,n_t); zeros(1,n_t)];
            end
        case 11
            % Seven-phase system, DOF: V0
            fi = 2/7*pi;
            A = 2/7*[1, cos(fi), cos(2*fi), cos(3*fi), cos(-3*fi), cos(-2*fi), cos(-fi);...
                0, sin(fi), sin(2*fi), sin(3*fi), sin(-3*fi), sin(-2*fi), sin(-fi);...
                1, cos(3*fi), cos(-fi), cos(2*fi), cos(-2*fi), cos(fi), cos(-3*fi);...
                0, sin(3*fi), sin(-fi), sin(2*fi), sin(-2*fi), sin(fi), sin(-3*fi);...
                1, cos(-2*fi), cos(3*fi), cos(fi), cos(-fi), cos(-3*fi), cos(2*fi);...
                0, sin(-2*fi), sin(3*fi), sin(fi), sin(-fi), sin(-3*fi), sin(2*fi)];
            if y_original
                ys = 1*[cos(wt); sin(wt); 0.15*cos(3*wt); 0.15*sin(3*wt); 0.08*cos(5*wt); 0.08*sin(5*wt)];
            else
                ys = 1*[cos(wt); sin(wt); zeros(1,n_t); zeros(1,n_t); zeros(1,n_t); zeros(1,n_t)];
            end
        case 12
            % Seven-phase system, fault-tolerant mode, DOF: Vx1, Vy1, Vx2, Vy2, V0
            fi = 2/7*pi;
            A = 2/7*[0, cos(fi), cos(2*fi), cos(3*fi), cos(-3*fi), cos(-2*fi), cos(-fi);...
                0, sin(fi), sin(2*fi), sin(3*fi), sin(-3*fi), sin(-2*fi), sin(-fi);];
            ys = 1*[cos(wt); sin(wt)];
        case 13
            % Nine-phase system, DOF: Vx1, Vy1, Vx2, Vy2, Vx3, Vy3, V0
            fi = 2/9*pi;
            A = 2/9*[1, cos(fi), cos(2*fi), cos(3*fi), cos(4*fi), cos(-4*fi), cos(-3*fi), cos(-2*fi), cos(-fi);...
                0, sin(fi), sin(2*fi), sin(3*fi), sin(4*fi), sin(-4*fi), sin(-3*fi), sin(-2*fi), sin(-fi)];
            ys = 1*[cos(wt); sin(wt)];
        case 14
            % Nine-phase system, DOF: Vx2, Vy2, Vx3, Vy3, V0
            fi = 2/9*pi;
            A = 2/9*[1, cos(fi), cos(2*fi), cos(3*fi), cos(4*fi), cos(-4*fi), cos(-3*fi), cos(-2*fi), cos(-fi);...
                0, sin(fi), sin(2*fi), sin(3*fi), sin(4*fi), sin(-4*fi), sin(-3*fi), sin(-2*fi), sin(-fi);...
                1, cos(3*fi), cos(-3*fi), 1, cos(3*fi), cos(-3*fi), 1, cos(3*fi), cos(-3*fi);...
                0, sin(3*fi), sin(-3*fi), 0, sin(3*fi), sin(-3*fi), 0, sin(3*fi), sin(-3*fi)];
            if y_original
                ys = 1*[cos(wt); sin(wt); 0.15*cos(3*wt); 0.15*sin(3*wt)];
            else
                ys = 1*[cos(wt); sin(wt); zeros(1,n_t); zeros(1,n_t)];
            end
        case 15
            % Nine-phase system, DOF: Vx3, Vy3, V0
            fi = 2/9*pi;
            A = 2/9*[1, cos(fi), cos(2*fi), cos(3*fi), cos(4*fi), cos(-4*fi), cos(-3*fi), cos(-2*fi), cos(-fi);...
                0, sin(fi), sin(2*fi), sin(3*fi), sin(4*fi), sin(-4*fi), sin(-3*fi), sin(-2*fi), sin(-fi);...
                1, cos(3*fi), cos(-3*fi), 1, cos(3*fi), cos(-3*fi), 1, cos(3*fi), cos(-3*fi);...
                0, sin(3*fi), sin(-3*fi), 0, sin(3*fi), sin(-3*fi), 0, sin(3*fi), sin(-3*fi);...
                1, cos(-4*fi), cos(fi), cos(-3*fi), cos(2*fi), cos(-2*fi), cos(3*fi), cos(-fi), cos(4*fi);...
                0, sin(-4*fi), sin(fi), sin(-3*fi), sin(2*fi), sin(-2*fi), sin(3*fi), sin(-fi), sin(4*fi);...
                ];
            if y_original
                ys = 1*[cos(wt); sin(wt); 0.15*cos(3*wt); 0.15*sin(3*wt); 0.08*cos(5*wt); 0.08*sin(5*wt)];
            else
                ys = 1*[cos(wt); sin(wt); zeros(1,n_t); zeros(1,n_t); zeros(1,n_t); zeros(1,n_t)];
            end
        case 16
            % Nine-phase system, DOF: V0
            fi = 2/9*pi;
            A = 2/9*[1, cos(fi), cos(2*fi), cos(3*fi), cos(4*fi), cos(-4*fi), cos(-3*fi), cos(-2*fi), cos(-fi);...
                0, sin(fi), sin(2*fi), sin(3*fi), sin(4*fi), sin(-4*fi), sin(-3*fi), sin(-2*fi), sin(-fi);...
                1, cos(3*fi), cos(-3*fi), 1, cos(3*fi), cos(-3*fi), 1, cos(3*fi), cos(-3*fi);...
                0, sin(3*fi), sin(-3*fi), 0, sin(3*fi), sin(-3*fi), 0, sin(3*fi), sin(-3*fi);...
                1, cos(-4*fi), cos(fi), cos(-3*fi), cos(2*fi), cos(-2*fi), cos(3*fi), cos(-fi), cos(4*fi);...
                0, sin(-4*fi), sin(fi), sin(-3*fi), sin(2*fi), sin(-2*fi), sin(3*fi), sin(-fi), sin(4*fi);...
                1, cos(-2*fi), cos(-4*fi), cos(3*fi), cos(fi), cos(-fi), cos(-3*fi), cos(4*fi), cos(2*fi);...
                0, sin(-2*fi), sin(-4*fi), sin(3*fi), sin(fi), sin(-fi), sin(-3*fi), sin(4*fi), sin(2*fi);...
                ];
            if y_original
                ys = 1*[cos(wt); sin(wt); 0.15*cos(3*wt); 0.15*sin(3*wt); 0.08*cos(5*wt); 0.08*sin(5*wt); 0.05*cos(7*wt); 0.05*sin(7*wt)];
            else
                ys = 1*[cos(wt); sin(wt); zeros(1,n_t); zeros(1,n_t); zeros(1,n_t); zeros(1,n_t); zeros(1,n_t); zeros(1,n_t)];
            end
        case 17
            % Nine-phase system, fault-tolerant mode, DOF: Vx1, Vy1, Vx2, Vy2, Vx3, Vy3, V0
            fi = 2/9*pi;
            A = 2/9*[0, cos(fi), cos(2*fi), cos(3*fi), cos(4*fi), cos(-4*fi), cos(-3*fi), cos(-2*fi), cos(-fi);...
                0, sin(fi), sin(2*fi), sin(3*fi), sin(4*fi), sin(-4*fi), sin(-3*fi), sin(-2*fi), sin(-fi)];
            ys = 1*[cos(wt); sin(wt)];
        case 18
            % Triple three-phase inverter for new compensator
            fi = 2/3*pi;
            A1 = 2/3*[1, cos(fi), cos(-fi);...
                0, sin(fi), sin(-fi);...
                1/2, 1/2, 1/2];
            A2 = [ 1  0 -1  0 -1  1  0  0  0;... % ua = u11 - u13 + u23 - u22
                0  0  0  1  0 -1  0 -1  1;... % ub = u21 - u23 + u33 - u32
                0 -1  1  0  0  0  1  0 -1;];  % uc = u31 - u33 + u13 - u12
            A = A1*A2;
            if y_original
                ys = 1*[cos(wt); sin(wt); -cos(wt)];
            else
                ys = 1*[cos(wt); sin(wt); zeros(1,n_t)];
            end
        otherwise
            disp("Error: System matrix A not defined.");
    end

    file_name = 'data/test_pars_' + string(j) + '_' + string(remove_same_columns) + '.mat';
    if isfile(file_name)
        load(file_name, "pars");
    else
        pars = Pars(A, remove_same_columns);
        save(file_name, "pars");
    end
    solver = Solver(pars);
    [n_y, n_x] = size(A);
    D = A' / (A*A');

    xs = zeros(n_x, n_t);
    xs_l2 = zeros(n_x, n_t);
    for k = 1:n_t
        xs(:,k) = solver.min_effort(ys(:,k));
        xs_l2(:,k) = D*ys(:,k);
    end

    ylims = [-max(xs_l2(:)); max(xs_l2(:))];

    fig = figure();
    set(gcf, 'Position', get(0, 'Screensize'));

    subplot(1, 3, 1);
    plot(ts, ys');
    xlabel('Time [s]');
    ylabel('Required voltage');
    title('Case ' + string(j))
    subplot(1, 3, 2);
    plot(ts, xs');
    xlabel('Time [s]');
    ylabel('Input voltage');
    title('Our approach');
    ylim(ylims);
    subplot(1, 3, 3);
    plot(ts, xs_l2');
    xlabel('Time [s]');
    ylabel('Input voltage');
    title('l2 approach');
    ylim(ylims);
end






