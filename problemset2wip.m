%% magic squares
T=zeros(1,100);
dvec=zeros(1,100);
T(1)=cputime;
for ii=1:100
    A=magic(100*ii);
    t1=cputime;
    dvec(ii)=det(A);
    T(ii)=cputime-t1;
end
%% plotting for problem 1
n=1:100;
n3=n.^3;
c=1.1*10^-4;
plot(c*n3,T);
hold on
plot(c*n3,c*n3);
hold off
grid on

%% Problem 3d

eps_vals = 10.^(-(0:9));
relerr_crout = zeros(size(eps_vals));
relerr_lu = zeros(size(eps_vals));
x_true = [1;1;1];

for k = 1:length(eps_vals)
    eps = eps_vals(k);
    A = [2 -2 0; eps-2 2 0; 0 -1 3];
    b = A*x_true;

    % Crout LU
    [L,U] = crout(A);
    y = L\b;
    x_crout = U\y;

    % matlab LU
    [L2,U2] = lu(A);
    x_lu = U2\(L2\b);

    relerr_crout(k) = norm(x_crout - x_true)/norm(x_true);
    relerr_lu(k) = norm(x_lu - x_true)/norm(x_true);
end

loglog(eps_vals, relerr_crout, 'o-', eps_vals, relerr_lu, 's-')
xlabel('\epsilon')
ylabel('Relative Error')
legend('Crout','MATLAB lu','Location','NorthWest')
title('Relative Error vs \epsilon for x = (1,1,1)^T')
grid on

%% Problem 3(e)

eps_vals = (1/3)*10.^(-(0:9));
relerr_crout = zeros(size(eps_vals));
relerr_lu = zeros(size(eps_vals));

x_true = [log(2.5);1;1];

for k = 1:length(eps_vals)
    eps = eps_vals(k);
    A = [2 -2 0; eps-2 2 0; 0 -1 3];
    b = A*x_true;

    % Crout LU
    [L,U] = crout(A);
    y = L\b;
    x_crout = U\y;

    % MATLAB LU
    [L2,U2] = lu(A);
    x_lu = U2\(L2\b);

    relerr_crout(k) = norm(x_crout - x_true)/norm(x_true);
    relerr_lu(k) = norm(x_lu - x_true)/norm(x_true);
end

loglog(eps_vals, relerr_crout, 'o-', eps_vals, relerr_lu, 's-')
xlabel('\epsilon')
ylabel('Relative Error')
legend('Crout','MATLAB lu','Location','NorthWest')
title('Relative Error vs \epsilon for x = (log(2.5),1,1)^T')
grid on
%% Problem 4
A=[1,1+0.5*10^-15,3;2,2,20;3,6,4];
[L,U]=crout(A);
R=A-L*U;
[L2 U2 P2]=lu(A);
RP=P2*A-L2*U2;

