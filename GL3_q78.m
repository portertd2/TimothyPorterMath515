% Question 7 and 8 script

% Example dataset: monthly average temperature (F)
x = 2006:2025; 
y = [0.64 0.67 0.54 0.66 ...
0.73 0.61 0.65 0.68 0.75 0.90 1.01 0.92 0.85 0.98 ...
1.01 0.85 0.90 1.17 1.28 1.19];
% create many points
x_dense = linspace(min(x),max(x),1450);

% all data
y_interp = cubic_spline(x,y,x_dense);

% every other point
x_half = x(1:2:end);
y_half = y(1:2:end);
y_half_interp = cubic_spline(x_half,y_half,x_dense);

% every 3rd point
x_third = x(1:3:end);
y_third = y(1:3:end);
y_third_interp = cubic_spline(x_third,y_third,x_dense);

%% Plot everything together

figure
hold on
% all data spline
plot(x_dense,y_interp,'b','LineWidth',1)

% every other point spline
plot(x_dense,y_half_interp,'r','LineWidth',1)

% every third point spline
plot(x_dense,y_third_interp,'g','LineWidth',1)

% original data points
plot(x,y,'ko-','MarkerSize',4,'LineWidth',1)

legend('Spline (all data)','Spline (every other point)','Spline (every 3rd point)','Original data','Location','NorthWest')

xlabel('Year')
ylabel('Global Land-Ocean Temperature Index (C)')
title('Cubic Spline Interpolation Comparison')

grid on