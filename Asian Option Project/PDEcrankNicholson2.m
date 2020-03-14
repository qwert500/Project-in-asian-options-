function [time, space, sol]=PDEcrankNicholson2(T,Z,n,m,r,sigma)

%======Initialization======
dt=T/n;
dz=2*Z/m;
d=dt/(dz^2);
sol=zeros(n+1,m+1);
B1=zeros(m+1,m+1,n+1);
B2=zeros(m+1,m+1,n+1);
time=zeros(1,n+1);
space=zeros(1,m+1);
space(1)=-Z;
A=zeros(m+1,m+1,n+1);

for i=2:n+1
  time(i)=time(i-1)+dt;
end

for j=2:m+1
  space(j)=space(j-1)+dz;
end

%=======Inital values================
for j=1:m+1
  sol(1,j)=max(space(j),0);
end
%=======boundary conditions===========
sol(:,1)=0;
sol(:,m+1)=Z;

%========Solving PDE==================
for i=1:n+1
  for k=2:m %diffrent for all t 
    A(k,k-1,i)=-d/2*sigma^2*((1-exp(-r*time(i)))/(r*T)-space(k))^2;
    A(k,k+1,i)=-d/2*sigma^2*((1-exp(-r*time(i)))/(r*T)-space(k))^2;
    A(k,k,i)=d*sigma^2*((1-exp(-r*time(i)))/(r*T)-space(k))^2;
  end
end

for i=1:n+1
  B1(:,:,i)=eye(m+1)+A(:,:,i)/2;
  B2(:,:,i)=eye(m+1)-A(:,:,i)/2;
end

B1(1,1,:)=1;
B1(m+1,m+1,:)=1;% For boundary conditions to be consistent
B2(1,1,:)=1;
B2(m+1,m+1,:)=1;

for i=2:n+1
  sol(i,:)=B1(:,:,i)\(B2(:,:,i-1)*sol(i-1,:)'); % Crank-Nicholson solution
end

end
