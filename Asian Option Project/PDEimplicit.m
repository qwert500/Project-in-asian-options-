function [time, space, sol]=PDEimplicit(T,Z,n,m,r,sigma)
%=====Initialization====
dt=T/n;
dz=Z/m;
d=dt/(dz^2);
sol=zeros(n+1,m+1);
time=zeros(1,n+1);
space=zeros(1,m+1);
A=zeros(m+1,m+1);

for i=2:n+1;
  time(i)=time(i-1)+dt;
end
space(1)=-Z;
for j=2:m+1
  space(j)=space(j-1)+2*dz;
end
%=====Determining A matrix=======
A(:,m+1)=Z;
A(:,1)=0;
%=====Initial solution===========
for j=1:m+1
sol(1,j)=max(space(j),0);
end
%=====Boundary conditions========
sol(:,1)=0;
sol(:,m+1)=Z;
%=====Solving PDE================
% gamma for t(i+1)
for i=2:n+1
  for k=2:m %diffrent for all t
    A(k,k-1)=-d/2*sigma^2*((1-exp(-r*time(i)))/(r*T)-space(k))^2;
    if k<m
    A(k,k+1)=-d/2*sigma^2*((1-exp(-r*time(i)))/(r*T)-space(k))^2;
    end
    A(k,k)=1+d*sigma^2*((1-exp(-r*time(i)))/(r*T)-space(k))^2;
  end
  A(1,1)=1;
  A(m+1,m+1)=1;
  sol(i,:)=sol(i-1,:)*transpose(inv(A));
end


