
function q=quantile(X,p)
% q=quantile(X,p)
% Compute pth quantile for each column in a matrix.
%
%  See Conover, "Practical Nonparametric Statistics", 2nd edition,p. 71
%
%***********  INPUT ARGUMENTS
%
% X (m1 x n1) matrix with variables as columns, observations as rows
% p (1 x 1)   the fraction of Xs smaller than q is <= p
%              the fraction of Xs greater than q is <= (1-p)
%
%
%*********  OUTPUT ARGUMENT
%
% q (1 x n1)  the values corresp to the pth quantile

Xs=sort(X);
[m1,n1]=size(X);


N=((1:m1)/(m1))';
Xss=[N Xs];
Xrev=N(m1:-1:1);

ilo = find(Xss(:,1)<=p);
ihi = find (Xrev <= (1-p));

if (length(ilo)+length(ihi)) == m1;
	q=(Xss(max(ilo),2:n1+1) + Xss(min(ihi),2:n1+1)) / 2.0;
else
	nout=m1-(length(ilo)+length(ihi));
	jlo=max(ilo)+1;
	jhi=max(ilo)+nout;
	q=(Xss(jlo,2:n1+1) + Xss(jhi,2:n1+1)) / 2.0;
end

if isempty(q)
	error('p must be  between 1/n and n/(n+1)')
end
