function [P,m] = pctn(X)

%  Transform a data array X into pctg of normal by dividing 
%  values by their col means.  Returtn the col means in m.

[m1,n1]=size(X);
m=mean(X);
M=m(ones(m1,1),:);
P=X ./ M;
