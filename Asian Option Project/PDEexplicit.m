function [time, space, sol]=PDEexplicit(T,Z,n,m,r,sigma)
%======Initialization======
dt=T/n;
dz=Z/m;
d=dt/(dz^2);
sol=zeros(n+1,m+1);
time=zeros(1,n+1);
space=zeros(1,m+1);
space(1)=-Z;
for i=2:n+1
  time(i)=time(i-1)+dt;
end

for j=2:m+1
  space(j)=space(j-1)+2*dz;
end
%=======Inital values================
for j=1:m+1
  sol(1,j)=max(space(j),0); 
end
%=======boundary conditions===========
sol(:,1)=0;
sol(:,m+1)=Z;


%========Solving PDE==================
%gamma at at i-1
for i=2:n+1
  for j=3:m+1
    sol(i,j-1)=(d*sigma^2)/2*(((1-exp(-r*time(i-1))))/(r*T)-space(j-1))^2*...
      (sol(i-1,j)-2*sol(i-1,j-1)+sol(i-1,j-2))+sol(i-1,j-1);
  end
end


  
  

