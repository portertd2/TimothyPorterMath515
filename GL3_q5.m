% Test script for Question 5
x = [0.005 0.010 0.020 0.050 0.100 0.200 0.500 1.000 2.000];
y = [0.924 0.896 0.859 0.794 0.732 0.656 0.536 0.430 0.316];
xstar = 0.032;
ystar = cubic_spline(x,y,xstar);
disp('Interpolated value:')
disp(ystar)