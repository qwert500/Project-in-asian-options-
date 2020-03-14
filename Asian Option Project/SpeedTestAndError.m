%=====Speed test=======

%=====Initialization======
clc, clear all, clear figures
T=1;
Z=5;
N=100000;
n=150;
m=1000; 
r=0.2;
d=(T/n)/((Z/m)^2);
sigma=1;
S0=100;
K=120;
%=====Finite difference method=======
tic
[time, space, sol_CN]=PDEcrankNicholson2(T,Z,n,m,r,sigma);
toc
%=====Monte carlo====================
tic
[price_MonteCarlo, conf95]=MonteCarlo_AC(S0,sigma,r,K,T,n,N);
toc
%=====ERROR=====
z0=1/(r*T)*(1-exp(-r*T))+(exp(-r*T))*(-K/S0);
[~,index]=min(abs(space-z0));
priceOfAsianCallOption_crankNicholson=S0*sol_CN(n+1,index);
error=100-100*(priceOfAsianCallOption_crankNicholson/price_MonteCarlo); % in procent



