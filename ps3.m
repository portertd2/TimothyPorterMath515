% MATH 515
% Problem Set 3 MATLAB Code
% Timothy Porter

%% Problem 3(a): Neville's Method
x = [0 1 2];
y = [2 -1 4];
xbar = 1.3;

p_neville = neville(x,y,xbar);
disp('Problem 3(a) Neville Approximation:')
disp(p_neville)

%% Problem 3(b): Newton Divided Differences
[dd,p_newton] = newtondd(x,y,xbar);

disp('Problem 3(b) Divided Difference Table:')
disp(dd)

disp('Problem 3(b) Newton Approximation:')
disp(p_newton)

%% Problem 4
f = @(x) x.*exp(-x);
a = -1;
b = 3;

% Uniform nodes
xU = linspace(a,b,5);

% Chebyshev nodes
k = 1:5;
xC = (a+b)/2 + (b-a)/2*cos((2*k-1)*pi/10);

% Legendre nodes (roots of P_5 scaled to [-1,3])
xLhat = [-0.9061798459 -0.5384693101 0 0.5384693101 0.9061798459];
xL = (a+b)/2 + (b-a)/2*xLhat;

% Interpolating polynomials
pU = polyfit(xU,f(xU),4);
pC = polyfit(xC,f(xC),4);
pL = polyfit(xL,f(xL),4);

% grid
xx = linspace(a,b,2000);

errU = f(xx) - polyval(pU,xx);
errC = f(xx) - polyval(pC,xx);
errL = f(xx) - polyval(pL,xx);

LinftyU = max(abs(errU));
LinftyC = max(abs(errC));
LinftyL = max(abs(errL));

L2U = sqrt(trapz(xx,errU.^2));
L2C = sqrt(trapz(xx,errC.^2));
L2L = sqrt(trapz(xx,errL.^2));

disp('Problem 4 Errors [Linfty, L2]:')
disp('Uniform:')
disp([LinftyU L2U])
disp('Chebyshev:')
disp([LinftyC L2C])
disp('Legendre:')
disp([LinftyL L2L])

%% Problem 5(a)
disp('Problem 5(a) Forward First Derivative Coefficients:')
disp(fdcoeff(1,[0 1]))

disp('Problem 5(a) Centered First Derivative Coefficients:')
disp(fdcoeff(1,[-1 1]))

disp('Problem 5(a) Centered Second Derivative Coefficients:')
disp(fdcoeff(2,[-1 0 1]))

%% Problem 5(b)
f = @(x) exp(x);
hvals = 10.^(-(1:6));
err = zeros(size(hvals));

for j = 1:length(hvals)
    h = hvals(j);
    approx = (f(h)-f(-h))/(2*h);
    err(j) = abs(approx - 1);
end

disp('Problem 5(b) Errors:')
disp([hvals' err'])

figure;
loglog(hvals,err,'o-')
xlabel('h')
ylabel('Error')
title('Problem 5(b): Centered Difference Convergence')
grid on
saveas(gcf,'problem5b.png')

%% Problem 5(c)
f = @(x) sin(x);
hvals = 10.^(-(0:5));
err = zeros(size(hvals));

for j = 1:length(hvals)
    h = hvals(j);
    dx = h*[0 1 2 3 4];
    c = fdcoeff(1,dx);
    approx = c' * f(dx)';
    err(j) = abs(approx - 1);
end

disp('Problem 5(c) Errors:')
disp([hvals' err'])

h = 1;
dx = h*[0 1 2 3 4];
c = fdcoeff(1,dx);
approx_h1 = c' * sin(dx)';
disp('Problem 5(c) Approximation for h=1:')
disp(approx_h1)

%% Problem 6(c)
f = @(x) log(x.^2 + 1);
x0 = 1;

h = 1;
D = zeros(4,4);

for i = 1:4
    D(i,1) = (f(x0) - f(x0-h))/h;
    h = h/2;
end

for j = 2:4
    for i = j:4
        D(i,j) = D(i,j-1) + ...
        (D(i,j-1)-D(i-1,j-1))/(2^(j-1)-1);
    end
end

disp('Problem 6(c) Richardson Table:')
disp(D)

%% Functions

function p = neville(x,y,xbar)
n = length(x);
Q = zeros(n,n);
Q(:,1) = y(:);

for j = 2:n
    for i = 1:n-j+1
        Q(i,j) = ((xbar - x(i))*Q(i+1,j-1) ...
        - (xbar - x(i+j-1))*Q(i,j-1)) ...
        / (x(i+j-1) - x(i));
    end
end

p = Q(1,n);
end

function [dd,pval] = newtondd(x,y,xbar)
n = length(x);
dd = zeros(n,n);
dd(:,1) = y(:);

for j = 2:n
    for i = 1:n-j+1
        dd(i,j) = (dd(i+1,j-1)-dd(i,j-1)) / (x(i+j-1)-x(i));
    end
end

pval = dd(1,1);
prodterm = 1;
for j = 2:n
    prodterm = prodterm * (xbar - x(j-1));
    pval = pval + dd(1,j) * prodterm;
    end
end

function c = fdcoeff(n,dx)
k = length(dx);
A = zeros(k,k);

for i = 1:k
    A(i,:) = dx.^(i-1);
end

b = zeros(k,1);
b(n+1) = factorial(n);

c = A\b;
end
