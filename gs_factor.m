function [Q,R] = gs_factor(A)
% [Q,R] =gs_factor(A)
% Modified Gram-Schmidt factorization.
% Produces A=QR, Q is mxn as is A
% R is nxn and upper triangular
% if R(k,k) is 0, q_k is set to 0 vector
% Timothy Porter, 2026-02-12

[m,n]=size(A);

Q=zeros(m,n);
R=zeros(n,n);

for k=1:n
    v=A(:,k);
    for l = 1:k-1
        R(l,k)=Q(:,l)'*v;
        v=v-R(l,k)*Q(:,l);
    end
    R(k,k)=norm(v);
    if R(k,k)>0
        Q(:,k)=v/R(k,k);
    else
        Q(:,k)=zeros(m,1);
    end
end

end
