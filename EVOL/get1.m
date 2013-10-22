function [B,yrs] = get1(X,c)

mvc = -999;
L1 = X(:,c+1)~=mvc;
t = X(L1,1);
x = X(L1,c+1);
tm=[min(t) max(t)];
n=length(x);
d = diff(tm)+1;
disp(['Series ',int2str(c)]);
disp(['Period: ',int2str(tm(1)),' ',int2str(tm(2))]);
disp(['Series Length ',int2str(n)]);
if d~=n
  error('Missing value internal to time series')
end

B=[t x];
yrs=tm;
