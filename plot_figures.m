function h = plot_figures(ts, xs, A1, A2, Um, y_lim, kappa)
    disc = linspace(0, 2*pi, 100);
    
    h = figure('visible', 'off', 'position', 4*[0, 0, 400, 100]);
    axis tight manual % this ensures that getframe() returns a consistent size
    
    subplot(1,4,1);
    plot(ts, xs');
    ylim(y_lim);
    ylabel('kappa = ' + string(kappa))
    title('xs');
    
    subplot(1,4,2);
    plot(ts, (A1*xs)');
    ylim(y_lim);
    title('A*xs');
    
    subplot(1,4,3);
    hold on;
    plot(ts, (A2([1 2],:)*xs)');
    plot(ts, Um*kappa*ones(1,length(ts)), 'k--');
    plot(ts, -Um*kappa*ones(1,length(ts)), 'k--');
    ylim(y_lim);
    title('B*xs');
    
    subplot(1,4,4);
    hold on;
    plot((A1(1,:)*xs)', (A1(2,:)*xs)');
    plot((A2(1,:)*xs)', (A2(2,:)*xs)');
    plot(Um*kappa*cos(disc), Um*kappa*sin(disc), 'k--');
    xlim(y_lim);
    ylim(y_lim);
    title('B*xs');
    
    drawnow;
end