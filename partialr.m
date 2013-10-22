function r=partialr(X1,X2);
% partialr: partial correlations for matrix of variables
% CALL: r=partialr(X1,X2);
%
% Meko 3-19-93
%
% Partial correlation coefficients
% Partials are computed between X1(:,1) and all the remaining cols of
% X1, removing the effects of all variables in X2
%
% Coded from equations in Mardia et al. 1979, p. 169-170

[m1,n1]=size(X1);
[m2,n2]=size(X2);
if m1~=m2;
	error('X1 AND X2 DO NOT HAVE SAME NUMBER OF ROWS!');
end

R=corrcoef([X1 X2]);  % correlation matrix between all n1+n2 series
[m3,n3]=size(R);

R22=R(n1+1:n3,n1+1:n3);
R1=R(n1+1:n3,1:n1);

rank(R22);
rstar=R(1:n1,1:n1) - R1' * inv(R22) * R1;

G=(diag(rstar)) * (diag(rstar))';
G1=G(1,:);
r= rstar(1,:) ./ sqrt(G1);


