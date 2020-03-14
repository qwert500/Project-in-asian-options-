function [time, space, sol]=PDEcrankNicholson(T,Z,n,m,r,sigma)

%======Initialization======
dt=T/n;
dz=Z/m;
d=dt/(dz^2);
sol=zeros(n+1,m+1);
solimp=zeros(n+1,m+1);
solexp=zeros(n+1,m+1);
time=zeros(1,n+1);
space=zeros(1,m+1);
space(1)=-Z;
A=zeros(m+1,m+1);
for i=2:n+1
  time(i)=time(i-1)+dt;
end

for j=2:m+1
  space(j)=space(j-1)+2*dz;
end

%=======Inital values================
for j=1:m+1
  sol(1,j)=max(space(j),0);
  solexp(1,j)=max(space(j),0);
  solimp(1,j)=max(space(j),0);
end
%=======boundary conditions===========
sol(:,1)=0;
solexp(:,1)=0;
solimp(:,1)=0;
sol(:,m+1)=Z;
solexp(:,m+1)=Z;
solimp(:,m+1)=Z;

%========Solving PDE==================

for i=2:n+1
  
  for j=3:m+1 % explicit solution
    solexp(i,j-1)=(d*sigma^2)/2*(((1-exp(-r*time(i-1))))/(r*T)-space(j-1))^2*...
      (sol(i-1,j)-2*sol(i-1,j-1)+sol(i-1,j-2))+sol(i-1,j-1);
  end
  
  for k=2:m %diffrent for all t
    A(k,k-1)=-d/2*sigma^2*((1-exp(-r*time(i+1)))/(r*T)-space(k))^2;
    if k<m
      A(k,k+1)=-d/2*sigma^2*((1-exp(-r*time(i+1)))/(r*T)-space(k))^2;
    end
    A(k,k)=1+d*sigma^2*((1-exp(-r*time(i+1)))/(r*T)-space(k))^2;
  end
  
  A(1,1)=1;
  A(m+1,m+1)=1;
  solimp(i,:)=sol(i-1,:)*transpose(inv(A)); % implicit solution
  
  for j=1:m+1
  sol(i,j)=(solexp(i-1,j).*solimp(i-1,j)); % Crank-Nicholson solution
  end
  
end

end













