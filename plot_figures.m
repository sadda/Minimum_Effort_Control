function h = plot_figures(ts, xs, A1, A2, Um, y_lim, kappa, visible)
    if nargin < 8
        visible = 'off';
    end
    disc = linspace(0, 2*pi, 100);
    
    h = figure('visible', visible, 'position', 4*[0, 0, 400, 100]);
    set(gcf, 'color', 'white')
    axis tight manual% this ensures that getframe() returns a consistent size
    
    pos_x = 0.04;
    width_x = 0.2;
    offset_x = 0.05;
    
    for i = 1:4
        subplot('Position', [pos_x+(i-1)*(width_x+offset_x) 0.12 width_x 0.82]);
        if i == 1
            plot(ts, xs');
            ylim(y_lim);
            xlabel('time [s]');
            ylabel('input voltage');
            title('x');
        elseif i == 2
            plot(ts, (A1*xs)');
            ylim(y_lim);
            xlabel('time [s]');
            ylabel('V_a, V_b');
            title('A*x (in abc)');
        elseif i == 3
            plot(ts, (A2([1 2],:)*xs)');
            hold on;
            plot(ts, Um*kappa*ones(1,length(ts)), 'k--');
            plot(ts, -Um*kappa*ones(1,length(ts)), 'k--');
            ylim(y_lim);
            xlabel('time [s]');
            ylabel('V_a, V_b');
            title('B*x (in abc)');
        elseif i == 4
            plot((A1(1,:)*xs)', (A1(2,:)*xs)');
            hold on;
            plot((A2(1,:)*xs)', (A2(2,:)*xs)');
            plot(Um*kappa*cos(disc), Um*kappa*sin(disc), 'k--');
            xlim(y_lim);
            ylim(y_lim);
            xlabel('V_\alpha');
            ylabel('V_\beta');
            title('A*x, B*x (in \alpha\beta)');
        end
        drawnow;
    end