function [L,U] = crout(A)
% implements Crout's Method for LU decomp on matrix A
[r,~]=size(A);
for ii=1:r
    L(ii,1)=A(ii,1);
    U(ii,ii)=1;
end
for jj=2:r
    U(1,jj)=A(1,jj)/L(1,1);
end
for ii=2:r
    for jj=2:ii
        L(ii,jj)=A(ii,jj)-L(ii,1:jj-1)*U(1:jj-1,jj);
    end
    for jj=ii+1:r
        U(ii,jj)=(A(ii,jj)-L(ii,1:ii-1)*U(1:ii-1,jj))/L(ii,ii);
    end
end

end

