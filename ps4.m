%% Problem Set 4 Code
%% Problem 4e
clear all;
clc;
f = @(t,y) y.*cos(t);
yexact = @(t) exp(sin(t));
h = pi/4;
t = 0:h:2*pi;
N = length(t)-1;

% Forward Euler
yf = zeros(size(t));
yf(1) = 1;
for n = 1:N
    yf(n+1) = yf(n) + h*f(t(n), yf(n));
end

% Backward Euler
yb = zeros(size(t));
yb(1) = 1;
for n = 1:N
    yb(n+1) = yb(n)/(1 - h*cos(t(n+1)));
end

% Exact solution
tfine = linspace(0,2*pi,1000);
yfine = yexact(tfine);
figure(1);
plot(tfine, yfine, 'k-', 'LineWidth', 1.5)
hold on
plot(t, yf, 'o-', 'LineWidth', 1.2)
plot(t, yb, 's-', 'LineWidth', 1.2)
xlabel('t')
ylabel('y')
legend('Exact', 'Forward Euler', 'Backward Euler', 'Location', 'Best')
title('Problem 4e: Exact and Euler Approximations')
grid on

%% Problem 5b
clear all;
clc;
y_exact = 0.5*(exp(1) - sin(1) - cos(1));
Nh=2;
h=1/Nh;
t=0:h:1;
y=zeros(size(t));
y(1)=0;
for n=1:Nh
    y(n+1)=y(n)+h*(sin(t(n)))+y(n);
end
y_approx=y(end);
error=abs(y_exact-y_approx);
fprintf('Forward approximation at t=1: %.10f\n', y_approx);
fprintf('Exactvalue at t=1: %.10f\n', y_exact)
fprintf('Error: %.10f\n', error)


%% Problem 5c
clear all; 
clc;
yexact = 0.5*(exp(1) - sin(1) - cos(1));
Nh_vals = 2.^(1:20);
h_vals = 1./Nh_vals;
err_fe = zeros(size(Nh_vals));
for k = 1:length(Nh_vals)
    Nh = Nh_vals(k);
    h = 1/Nh;
    t = linspace(0,1,Nh+1);
    y = zeros(1,Nh+1);
    y(1) = 0;
    for n = 1:Nh
        y(n+1) = y(n) + h*(sin(t(n)) + y(n));
    end
    err_fe(k) = abs(yexact - y(end));
end
figure(2);
loglog(h_vals, err_fe, 'o-','LineWidth',1.2)
xlabel('h')
ylabel('Error at t = 1')
title('Problem 5c: Forward Euler Error')
grid on

%% Problem 5d
clear all; 
clc;
yexact = 0.5*(exp(1) - sin(1) - cos(1));
Nh_vals = 2.^(1:10);
h_vals = 1./Nh_vals;
err_be = zeros(size(Nh_vals));
for k = 1:length(Nh_vals)
    Nh = Nh_vals(k);
    h = 1/Nh;
    t = linspace(0,1,Nh+1);
    y = zeros(1,Nh+1);
    y(1) = 0;
    for n = 1:Nh
        y(n+1) = (y(n) + h*sin(t(n+1)))/(1-h);  
    end
    err_be(k) = abs(yexact - y(end));
end
figure(3);
loglog(h_vals, err_be, 's-','LineWidth',1.2)
xlabel('h')
ylabel('Error a t = 1')
title('Problem 5d: Backward Euler Error')
grid on