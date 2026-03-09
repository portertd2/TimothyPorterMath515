function ystar=cubic_spline(x,y,xstar)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Please copy this code into a new file and rename cublic_spline.m
% Change the function name in line 1 also to cubic_spline.

% Questions 1-3 ask you to sequentially add to the template to create a
% cubic spline interpolator. Much of the info is given to you but you will
% need to think beyond "cut and paste".
% Question 4 asks you to interpet code already given.
% Questions 5-6 ask you to run your interpolator with two different sets of
% data. You may want to create a script that includes both data sets, from
% which you run cubic_spline.m.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Question 1

% n relative to # of data points
n=length(x)-1;

%Initialize A and rvec (the RHS)
A=zeros(n+1,n+1);
rvec=zeros(n+1,1);

%Calculate vector of h values (not necessarily equally spaced)
h=zeros(1,n);
for jj=1:n
    h(jj)=x(jj+1)-x(jj);
end

%Populate upper left and lower right of A
A(1,1)=1;
A(n+1,n+1)=1;

%Populate rest of matrix A and answer vector r
for jj=2:n
    A(jj,jj-1)=h(jj-1);
    A(jj,jj)=2*(h(jj-1)+h(jj));
    A(jj,jj+1)=h(jj);
end
for ii=2:n
    rvec(ii)=3*(((y(ii+1)-y(ii))/h(ii))-((y(ii)-y(ii-1))/h(ii-1)));
end

%%% Question 2

c=A\rvec;

%%% Question 3
b=zeros(1,n);
d=zeros(1,n);

%(a)
a=y;

%(b-c)
for j=1:n
    d(j)=(c(j+1)-c(j))/(3*h(j));
    b(j)=(a(j+1)-a(j))/(h(j)) - (1/3)*(2*c(j)+c(j+1))*h(j);
end

%%% At this point you have calculated all values for a,b,c,d

%%% Question 4
% You do not need to add code -- the question asks you to determine
% what is happening in the existing code given to you here.

coeffs=[];
for k=1:n
    coeffs(k,:)=[d(k) c(k) b(k) a(k)];
end

breaks=[x];
pp=mkpp(breaks,coeffs);
ystar=ppval(pp,xstar);

%%% Question 5-6
% You do not need to add code here -- the questions ask you to interpolate
% two different sets of data and evaluate at two xstar values, one with a known output 
% for ystar to check your code. 