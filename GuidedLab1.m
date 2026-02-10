% Heron's method

promptA  = 'Input the number a to square root: ';
promptY0 = 'Input the initial guess y0: ';
promptN  = 'Input the number of iterations N: ';

a  = input(promptA);
y0 = input(promptY0);
N  = input(promptN);

% Exact solution (used as "true" value)
y_exact = sqrt(a);

% Preallocate vectors
yvec        = zeros(N+1,1);
true_relerr = zeros(N+1,1);
approx_relerr = zeros(N+1,1);

% Initial condition
yvec(1)        = y0;
true_relerr(1) = abs(y0 - y_exact)/y_exact;
approx_relerr(1) = NaN;   % not defined for n = 0

y = y0;

% Newton iteration
for ii = 2:N+1
    y_new = 0.5*(y + a/y);
    
    yvec(ii) = y_new;
    
    % True relative error
    true_relerr(ii) = abs(y_new - y_exact)/y_exact;
    
    % Approximate relative error
    approx_relerr(ii) = abs(y_new - y)/abs(y_new);
    
    y = y_new;
end

% Iteration index
ivec = 0:N;

% Log-log plot
figure
loglog(ivec, true_relerr, 'o-','LineWidth',1.5)
hold on
loglog(ivec(2:end), approx_relerr(2:end), 's--','LineWidth',1.5)
hold off

xlabel('Iteration')
ylabel('Relative Error')
legend('True relative error', 'Approximate relative error', ...
       'Location','southwest')
title('True vs Approximate Relative Error for Heron''s Method')
grid on

%% Method Comparison: Relative Error vs Iteration (log-log)

promptA  = 'Input the number a to square root: ';
promptY0 = 'Input the initial guess y0: ';
promptN  = 'Input the number of iterations N: ';

a  = input(promptA);
y0 = input(promptY0);
N  = input(promptN);

% Exact solution
y_exact = sqrt(a);

% Preallocate vectors
m1vec = zeros(N+1,1);  % Newton / Heron
m2vec = zeros(N+1,1);  % Unstable method
m3vec = zeros(N+1,1);  % Linear method
m4vec = zeros(N+1,1);  % Cubic method
m5vec = zeros(N+1,1);  % Additional method

m1err = zeros(N+1,1);
m2err = zeros(N+1,1);
m3err = zeros(N+1,1);
m4err = zeros(N+1,1);
m5err = zeros(N+1,1);

% Initial conditions
m1vec(1) = y0;
m2vec(1) = y0;
m3vec(1) = y0;
m4vec(1) = y0;
m5vec(1) = y0;

m1err(1) = abs(y0 - y_exact)/y_exact;
m2err(1) = abs(y0 - y_exact)/y_exact;
m3err(1) = abs(y0 - y_exact)/y_exact;
m4err(1) = abs(y0 - y_exact)/y_exact;
m5err(1) = abs(y0 - y_exact)/y_exact;

m1 = y0;
m2 = y0;
m3 = y0;
m4 = y0;
m5 = y0;

% Iterations
for ii = 2:N+1
    m1new = 0.5*(m1 + a/m1);                        % Newton (quadratic)
    m2new = m2 + m2^2 - a;                          % Unstable method
    m3new = 0.25*(3*m3 + 2/m3);                     % Linear method
    m4new = (1/8)*(3*m4 + 12/m4 - 4/m4^3);          % Cubic method
    m5new = m5*((m5^2 + 3*a)/(3*m5^2 + a));         % Additional method (Halley's)
    
    m1vec(ii) = m1new;
    m2vec(ii) = m2new;
    m3vec(ii) = m3new;
    m4vec(ii) = m4new;
    m5vec(ii) = m5new;
    
    m1err(ii) = abs(m1new - y_exact)/y_exact;
    m2err(ii) = abs(m2new - y_exact)/y_exact;
    m3err(ii) = abs(m3new - y_exact)/y_exact;
    m4err(ii) = abs(m4new - y_exact)/y_exact;
    m5err(ii) = abs(m5new - y_exact)/y_exact;
    
    m1 = m1new;
    m2 = m2new;
    m3 = m3new;
    m4 = m4new;
    m5 = m5new;
end

% Iteration index
ivec = 0:N;

% Log-log plot
figure
loglog(ivec, m1err, 'o-','LineWidth',1.5); hold on
loglog(ivec, m2err, 's--','LineWidth',1.5);
loglog(ivec, m3err, 'd-','LineWidth',1.5);
loglog(ivec, m4err, '^-','LineWidth',1.5);
loglog(ivec, m5err, 'x-','LineWidth',1.5);
hold off

xlabel('Iteration')
ylabel('Relative Error')
title(sprintf('Relative Error vs Iteration (a = %.4g, y_0 = %.4g, N = %d)', a, y0, N))
legend('Newton (quadratic)', ...
       'Unstable method', ...
       'Linear method', ...
       'Cubic method', ...
       'Halley''s method', ...
       'Location','northwest')
grid on

