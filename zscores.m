function z=zscores(x);
% z=zscores(x);
% zscores.m   convert a matrix to zscores.
% Feb 28, 1992, by D.Meko
%** INPUT ARGUMENT
%  x (? x ?)  input data, rows observations, 
%		columns variables. May alternatively be a rv.
%****  OUTPUT ARG
%
% z  (? x ?)  corresponding standardized variables -mean 0, st dev 1


[m1,n1]=size(x);

xm = mean(x);
xs = std (x);

if m1==1;  % Handle row-vector case

	m1=n1;  % number of years is number of elements
	z= (x - xm(:,ones(m1,1))) ./  xs(:,ones(m1,1));
else
	z= (x - xm(ones(m1,1),:)) ./  xs(ones(m1,1),:);
end
	
