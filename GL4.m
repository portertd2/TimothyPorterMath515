% MATH 515 - Guided Lab 4

clear; clc; close all;

function [t,w] = rk4_full(f,a,b,N,y0)
h = (b-a)/N;
t = linspace(a,b,N+1);
w = zeros(length(y0), N+1);
w(:,1) = y0;

for i = 1:N
    k1 = h*f(t(i), w(:,i));
    k2 = h*f(t(i)+h/2, w(:,i)+k1/2);
    k3 = h*f(t(i)+h/2, w(:,i)+k2/2);
    k4 = h*f(t(i)+h, w(:,i)+k3);

    w(:,i+1) = w(:,i) + (k1 + 2*k2 + 2*k3 + k4)/6;
end
end

function w_final = rk4_final(f,a,b,h,y0)
N = round((b-a)/h);
t = a;
w = y0;

for i = 1:N
    k1 = h*f(t,w);
    k2 = h*f(t+h/2, w+k1/2);
    k3 = h*f(t+h/2, w+k2/2);
    k4 = h*f(t+h, w+k3);

    w = w + (k1 + 2*k2 + 2*k3 + k4)/6;
    t = t + h;
end

w_final = w;
end

f = @(t,y) 1 - y + exp(2*t)*y.^2;
y_exact = @(t) exp(-t).*tan(exp(t)-1);

h = 0.1;
N = round(0.9/h);
[t,w] = rk4_full(f,0,0.9,N,0);

err = y_exact(t) - w;

figure;
plot(t, w, 'o-', t, y_exact(t), '-');
legend('RK4','Exact');
title('RK4 Approximation vs Exact Solution');

hs = logspace(-5,-1,20);

errors = zeros(size(hs));
h_actual = zeros(size(hs));

for k = 1:length(hs)

    h = hs(k);

    N = ceil(0.9/h);

    h_actual(k) = 0.9/N;

    w_final = rk4_final(f,0,0.9,h_actual(k),0);

    errors(k) = abs(y_exact(0.9) - w_final);

end

figure;
loglog(h_actual, errors, 'o-');
xlabel('h');
ylabel('Error');
title('RK4 Convergence Study');

% Use only h >= 1e-4 for slope estimation
idx = h_actual >= 1e-4;

p = polyfit(log(h_actual(idx)), log(errors(idx)), 1);

disp(['Estimated order (using h >= 1e-4): ', num2str(p(1))]);

disp(['Estimated order: ', num2str(p(1))]);

m = 16; a = 0.25;

f_sys = @(t,y) [
    1 - y(1) - (m*y(1)*y(2))/(a+y(1));
    (m*y(1)*y(2))/(a+y(1)) - y(2)
];

h = 0.05;
N = round(12/h);
[t,w] = rk4_full(f_sys,0,12,N,[0.5;0.02]);

figure;
plot(t,w(1,:), t,w(2,:));
legend('s(t)','x(t)');
title('Chemostat Model (h = 0.05)');

h = 0.06;
N = round(12/h);
[t,w] = rk4_full(f_sys,0,12,N,[0.5;0.02]);

figure;
plot(t,w(1,:), t,w(2,:));
legend('s(t)','x(t)');
title('Chemostat Model (h = 0.06)');

x = linspace(-3,3,300);
y = linspace(-3,3,300);
[X,Y] = meshgrid(x,y);
Z = X + 1i*Y;

Q = 1 + Z + Z.^2/2 + Z.^3/6 + Z.^4/24;

figure;
contour(X,Y,abs(Q),[1 1]);
xlabel('Re(z)'); ylabel('Im(z)');
title('|Q(z)| = 1 Stability Boundary');

lambda = -55.31;
f_test = @(t,y) lambda*y;

h1 = 0.01;
N1 = round(10/h1);
[t1,w1] = rk4_full(f_test,0,10,N1,1);

h2 = 0.1;
N2 = round(10/h2);
[t2,w2] = rk4_full(f_test,0,10,N2,1);

figure;
plot(t1,w1,'b', t2,w2,'r--');
legend('Stable h','Unstable h');
title('Stability Comparison');

lambda = 1 + 20i;
h = 0.1;
z = h*lambda;

figure;
contour(X,Y,abs(Q),[1 1]); hold on;
plot(real(z), imag(z), 'ro','MarkerSize',8,'LineWidth',2);
xlabel('Re(z)'); ylabel('Im(z)');
title('Stability Region with h\lambda Point');

f_test = @(t,y) lambda*y;
[t,w] = rk4_full(f_test,0,10,100,1);

figure;
plot(t, real(w));
title('Real Part of RK4 Solution (lambda = 1 + 20i)');