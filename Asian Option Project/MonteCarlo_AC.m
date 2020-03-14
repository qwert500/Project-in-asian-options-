function [price, conf95]=MonteCarlo_AC(s,sigma,r,K,T,N,n)

stockPath=StockPath(s,sigma,r,T,N,n);
payOff=max(0,mean(stockPath)-K);
price=exp(-r*T)*mean(payOff);
conf95=1.96*std(payOff)/sqrt(n);
