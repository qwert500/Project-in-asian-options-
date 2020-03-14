%=====Asian Option Project=====

%=====Input parameters=========
clc, clear all, clear figures
T=1;
Z=1;
n=100;
m=100; 
r=0.1;
d=(T/n)/((Z/m)^2)
sigma=0.5;

%=====Solving PDE==============
[time, space, sol_CN,B1,B2]=PDEcrankNicholson2(T,Z,n,m,r,sigma);


%=====Figures for crank-Nicholson method======
anim(space,sol_CN,0.8)
figure(2)
surf(time, space, transpose(sol_CN))
zlabel({'u(t,z)'})
xlabel({'t'})
ylabel({'z'})


