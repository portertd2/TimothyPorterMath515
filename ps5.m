%% Problem 2: Convergence comparison for second-order methods
clear; clc; close all;

% IVP:
% y' = 4t sqrt(1+y^2)/y, y(0)=1, 0 <= t <= 5

f = @(t,y) 4*t.*sqrt(1+y.^2)./y;
yexact = @(t) sqrt((2*t.^2 + sqrt(2)).^2 - 1);

a = 0;
b = 5;
alpha = 1;

% Step sizes
Nvals = [50 100 200 400 800 1600];
hvals = (b-a)./Nvals;

err_rk2 = zeros(size(hvals));
err_ab2 = zeros(size(hvals));
err_am1 = zeros(size(hvals));
err_heun = zeros(size(hvals));

for m = 1:length(Nvals)

    N = Nvals(m);
    h = hvals(m);
    t = linspace(a,b,N+1);

    %% (a) Optimal RK2 method
    w = zeros(1,N+1);
    w(1) = alpha;

    for j = 1:N
        k1 = f(t(j),w(j));
        k2 = f(t(j) + (2/3)*h, w(j) + (2/3)*h*k1);
        w(j+1) = w(j) + h*((1/4)*k1 + (3/4)*k2);
    end

    err_rk2(m) = max(abs(w - yexact(t)));

    %% (b) Two-step Adams-Bashforth method
    w = zeros(1,N+1);
    w(1) = alpha;

    % Use exact starting value at t = h
    w(2) = yexact(t(2));

    for j = 2:N
        w(j+1) = w(j) + (h/2)*(3*f(t(j),w(j)) - f(t(j-1),w(j-1)));
    end

    err_ab2(m) = max(abs(w - yexact(t)));

    %% (c) One-step Adams-Moulton method
    w = zeros(1,N+1);
    w(1) = alpha;

    for j = 1:N
        fj = f(t(j),w(j));

        % Solve implicit equation:
        % W = w_j + h/2 [f(t_j,w_j) + f(t_{j+1},W)]
        g = @(W) W - w(j) - (h/2)*(fj + f(t(j+1),W));

        % Predictor used as initial guess
        Wguess = w(j) + h*fj;

        w(j+1) = fzero(g,Wguess);
    end

    err_am1(m) = max(abs(w - yexact(t)));

    %% (d) Heun method: predictor-corrector
    w = zeros(1,N+1);
    w(1) = alpha;

    for j = 1:N
        k1 = f(t(j),w(j));

        % Predictor
        p = w(j) + h*k1;

        % Corrector
        k2 = f(t(j+1),p);
        w(j+1) = w(j) + (h/2)*(k1 + k2);
    end

    err_heun(m) = max(abs(w - yexact(t)));

end

%% Compute observed convergence rates
slope_rk2 = polyfit(log(hvals),log(err_rk2),1);
slope_ab2 = polyfit(log(hvals),log(err_ab2),1);
slope_am1 = polyfit(log(hvals),log(err_am1),1);
slope_heun = polyfit(log(hvals),log(err_heun),1);

fprintf('Observed slope for Optimal RK2: %.4f\n', slope_rk2(1));
fprintf('Observed slope for AB2: %.4f\n', slope_ab2(1));
fprintf('Observed slope for Adams-Moulton: %.4f\n', slope_am1(1));
fprintf('Observed slope for Heun: %.4f\n', slope_heun(1));

%% Plot error vs step size
figure;
loglog(hvals,err_rk2,'o-','LineWidth',1.5); hold on;
loglog(hvals,err_ab2,'s-','LineWidth',1.5);
loglog(hvals,err_am1,'^-','LineWidth',1.5);
loglog(hvals,err_heun,'d-','LineWidth',1.5);

grid on;
xlabel('Step size h');
ylabel('Maximum global error');
title('Problem 2: Error vs. Step Size');
legend('Optimal RK2','Two-step AB','One-step AM','Heun','Location','best');

%% Problem 4: Lorenz equations using ODE45
clear; clc; close all;

sigma = 10;
r = 28;
b = 8/3;

lorenz = @(t,u) [
    sigma*(u(2)-u(1));
    r*u(1) - u(2) - u(1)*u(3);
    u(1)*u(2) - b*u(3)
];

t0 = 0;
tf = 10;
h = 0.001;
tspan = t0:h:tf;

u0 = [5; 5; 5];
u0tilde = [5; 5; 5] + 1e-5*[1; 1; 1];

opts = odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',h);

%% Part (a): Original initial condition using ode45
[t,u] = ode45(lorenz,tspan,u0,opts);

x = u(:,1);
y = u(:,2);
z = u(:,3);

figure;
plot(x,z,'LineWidth',1);
grid on;
xlabel('x');
ylabel('z');
title('Lorenz System: x vs. z, Initial Condition (5,5,5)');

figure;
plot3(x,y,z,'LineWidth',1);
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
title('Lorenz System: 3D Plot, Initial Condition (5,5,5)');

%% Part (b): Perturbed initial condition using ode45
[ttilde,utilde] = ode45(lorenz,tspan,u0tilde,opts);

xtilde = utilde(:,1);
ytilde = utilde(:,2);
ztilde = utilde(:,3);

figure;
plot(xtilde,ztilde,'LineWidth',1);
grid on;
xlabel('x');
ylabel('z');
title('Lorenz System: x vs. z, Perturbed Initial Condition');

figure;
plot3(xtilde,ytilde,ztilde,'LineWidth',1);
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
title('Lorenz System: 3D Plot, Perturbed Initial Condition');


%% Part (c): L2 norm of difference between solutions
difference = u - utilde;
diffnorm = sqrt(sum(difference.^2,2));

figure;
semilogy(t,diffnorm,'LineWidth',1.5);
grid on;
xlabel('t');
ylabel('L^2 norm of solution difference');
title('Difference Between Lorenz Solutions with Nearby Initial Conditions');

%% Problem 4(d): Different solvers and step sizes
clear; clc; close all;
sigma = 10;
r = 28;
b = 8/3;
lorenz = @(t,u) [
    sigma*(u(2)-u(1));
    r*u(1) - u(2) - u(1)*u(3);
    u(1)*u(2) - b*u(3)
];
t0 = 0;
tf = 10;
u0 = [5; 5; 5];
u0tilde = [5; 5; 5] + 1e-5*[1; 1; 1];
solvers = {@ode45, @ode23, @ode113};
solver_names = {'ode45','ode23','ode113'};

hvals = [0.2, 0.1, 0.05, 0.025];
RelTol_value = 1e-4;
AbsTol_value = 1e-6;

num_runs = length(solvers)*length(hvals);
results_solver = strings(num_runs,1);
results_h = zeros(num_runs,1);
results_d10 = zeros(num_runs,1);
results_maxd = zeros(num_runs,1);
run_index = 1;

figure;

fprintf('\nProblem 4(d): Solver and step-size comparison\n');
fprintf('%8s %12s %15s %15s\n','Solver','MaxStep','d(10)','max d(t)');
fprintf('%8s %12s %15s %15s\n','------','-------','-----','--------');

for s = 1:length(solvers)
    solver = solvers{s};
    solver_name = solver_names{s};
    for k = 1:length(hvals)
        h = hvals(k);
        tspan = t0:h:tf;
        opts = odeset( ...
            'RelTol', RelTol_value, ...
            'AbsTol', AbsTol_value, ...
            'MaxStep', h);
        [tA,uA] = solver(lorenz,tspan,u0,opts);
        [tB,uB] = solver(lorenz,tspan,u0tilde,opts);
        diff = uA - uB;
        diffnorm = sqrt(sum(diff.^2,2));
        d10 = diffnorm(end);
        maxd = max(diffnorm);
        fprintf('%8s %12.4f %15.6e %15.6e\n', ...
            solver_name, h, d10, maxd);
        results_solver(run_index) = solver_name;
        results_h(run_index) = h;
        results_d10(run_index) = d10;
        results_maxd(run_index) = maxd;

        subplot(length(solvers),length(hvals),run_index);
        semilogy(tA,diffnorm,'LineWidth',1.2);
        grid on;
        xlabel('t');
        ylabel('d(t)');
        title(sprintf('%s, h = %.3f',solver_name,h));

        run_index = run_index + 1;

    end
end

sgtitle('Problem 4(d): Difference Norms for Different Solvers and Step Sizes');

results_table = table(results_solver, results_h, results_d10, results_maxd, ...
    'VariableNames', {'Solver','MaxStep','d_at_10','max_d'});

disp(' ');
disp('Summary table:');
disp(results_table);

