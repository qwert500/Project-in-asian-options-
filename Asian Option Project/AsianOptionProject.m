%=====Asian Option Project=====
tic
%=====Input parameters=========
clc, clear all, clear figures

T=1/12;
Z=5;
n=100;
N=300000; % # of geometric brownian motions for the monte carlo method...
m=100; 
r=0.01;
d=(T/n)/((Z/m)^2);
sigma=0;
K=90;
S0=100; % initial stock price
qmax=40;
sigmaStep=0.5;
priceOfEuropeanCall=zeros(qmax,1);
price_MonteCarlo=zeros(qmax,1);
priceOfAsianCallOption_crankNicholson=zeros(qmax,1);

Kmax=100;
putCallParity=zeros(Kmax,1);

for q=1:qmax
  sigma=sigma+sigmaStep;

%=====Solving PDE==============
[time, space, sol_CN]=PDEcrankNicholson2(T,Z,n,m,r,sigma);


%=====Figures for crank-Nicholson method======

%anim(space,sol_CN,0.8)
%figure(2)
%surf(time, space, transpose(sol_CN))
%zlabel({'u(t,z)'})
%xlabel({'t'})
%ylabel({'z'})

%=====Pricing the Asian Option with solution from Crank-Nicholson Method==

z0=1/(r*T)*(1-exp(-r*T))+(exp(-r*T))*(-K/S0);
[~,index]=min(abs(space-z0));

priceOfAsianCallOption_crankNicholson(q)=S0*sol_CN(n+1,index);

%====Pricing the the normal option =========
[priceOfEuropeanCall(q),~]=blsprice(S0,K,r,T,sigma);

%====Pricing the asian option with the monte carlo method=====
%[price_MonteCarlo(q), conf95]=MonteCarlo_AC(S0,sigma,r,K,T,n,N);

end
%====Price plot for the methods on AO and the standard price @ variance
figure(1)
hold on
plot(sigmaStep:sigmaStep:qmax*sigmaStep,priceOfEuropeanCall,'y',...
  sigmaStep:sigmaStep:qmax*sigmaStep,priceOfAsianCallOption_crankNicholson,'r')
xlabel('variance')
ylabel('Price of option')
title('Comparison of standard call options and asian call options')
legend('european call', 'Asian call, finite difference method')
toc
%%
clc, clear all 
%========Put-call pairity=========
% section for calculating the put-call parity by changing the loop variable

figure(2)
q=0;
S0=100;
K=90;
r=0.01;
T=1/12;
for r=0:0.01:0.4
  q=q+1;
putCallParity(q)=S0*(1-exp(-r*T))/(r*T)-K*exp(-r*T);
end
plot(0:0.01:0.4,putCallParity)
xlabel('r')
ylabel('put call pairity')
title('putcall pairity for the asian option')
%===================================


