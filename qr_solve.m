function x = qr_solve(Q, R, b)
% x = qr_solve(Q, R, b)
% solves the linear system Ax = b using a QR factorization A = QR
% Assumes: 
%  - Q is an orthogonal matrix (Q' * Q = I)
%  - R is an upper triangular matrix
%  - Q and R come from gs_factor or built-in qr
%2
% The method uses
%   Q R x = b
%   R x = Q' b
% and then solves the upper triangular system by back substitution.
%
% Timothy Porter
% 2026-02-12

y = Q'*b;

n = length(y);
x = zeros(n,1);

for i = n:-1:1
    x(i) = (y(i) - R(i,i+1:n)*x(i+1:n)) / R(i,i);
end

end
