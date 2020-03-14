function Path=StockPath(s,sigma,r,T,N,n)
h=T/N;
W=randn(n,N);
q=ones(n,N);
Path=s*exp((r-sigma^2/2)*h.*cumsum(q')+sigma*sqrt(h)*cumsum(W'));
Path=[s*ones(1,n);Path];