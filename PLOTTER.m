function [] = PLOTTER(x, xhat, attractor)

 if ~exist('attractor','var')
    attractor = [0;0]';
 end
if exist('xhat','var')
    figure, hold on, scatter(x(1,1),x(2,1), 'filled'),scatter(attractor(1),attractor(2), 'filled'), plot(x(1,:),x(2,:)), plot(xhat(1,:),xhat(2,:)), legend('$x_0$', 'Attractor', 'Trajectory', 'Position Estimate', location='best'), title('XY Coordinates');
    
else
    figure, hold on, scatter(x(1,1),x(2,1), 'filled'),scatter(attractor(1),attractor(2), 'filled'), plot(x(1,:),x(2,:)), legend('$x_0$', 'Attractor', 'Trajectory', location='best'), title('XY Coordinates');
end
figure, plot(x(3,:)), title('Robot Angle \theta');
end

