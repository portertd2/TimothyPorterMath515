function Q = gram_schmidt(X)
% Q = gram.schmidt(X)
% input the matrix X to obtain the matrix Q, the result of the Gram-Schmidt
% process, which generates a collection of orthonormal vectors spanning the
% same space as X
% Timothy Porter, 2026-02-12

[m,n]=size(X);
Q=zeros(m,n);
nq=0;
R=zeros(n,n);

for k=1:n
    for l = 1:nq
        R(l,k)=Q(:,l)'*X(:,k);
        X(:,k)=X(:,k)-R(l,k)*Q(:,l);
    end
    R(k,k)=norm(X(:,k));
    if R(k,k)>0
        nq=nq+1;
        Q(:,nq)=X(:,k)/R(k,k);
    end
end

Q =Q(:,1:nq);

end
