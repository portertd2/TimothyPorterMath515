%% Problem 1b
eps = 5e-3;
sum_pi = 0;
n = 0;

while 4/(2*n+3) > eps
    sum_pi = sum_pi + (-1)^n/(2*n+1);
    n = n + 1;
end

pi_approx = 4*sum_pi;
error = abs(pi - pi_approx);

fprintf('Approximation: %.6f\n', pi_approx)
fprintf('Actual error: %.6e\n', error)
fprintf('Number of terms: %d\n', n)
%% 1c
eps_vals = logspace(-4,-1,30);   % tolerances
N_vals = zeros(size(eps_vals));
err_act = zeros(size(eps_vals));

for k = 1:length(eps_vals)
    eps = eps_vals(k);
    sum_pi = 0;
    n = 0;

    while 4/(2*n+3) > eps
        sum_pi = sum_pi + (-1)^n/(2*n+1);
        n = n + 1;
    end

    pi_approx = 4*sum_pi;
    N_vals(k) = n;
    err_act(k) = abs(pi - pi_approx);
end

loglog(N_vals, eps_vals, 'ro-', ...
       N_vals, err_act, 'bx-')
grid on
xlabel('Number of terms N')
ylabel('Error / Tolerance')
legend('\epsilon (tolerance)', 'Actual error', 'Location','southwest')
%% 2
x = logspace(-6,-1,100);
f1 = sin(x.^2)./x.^2;
f2 = (sin(x)).^2./x.^2;

loglog(x, abs(f1-1), x, abs(f2-1))
legend('sin(x^2)/x^2','(sin x)^2/x^2')
xlabel('x')
ylabel('Error')
grid on
%% 5a
x = linspace(-5e-6,5e-6,1001);
f = log(1+x) - cos(x) - x + 1;
plot(x,f), grid on
title('Problem 5a');

%% 5b
f2 = x.^3/3 - x.^4/24;
plot(x,f2), grid on
title('Problem 5b');
%% 7
A = [3.4 2.8; 8 6.6];
b1 = [3.4; 8];
b2 = [3.41; 8];

x1 = A\b1
x2 = A\b2