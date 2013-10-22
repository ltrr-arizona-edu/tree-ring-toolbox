function [E,K1,W]=whit2(Y,nhi)

% Sept 20, 1991, by Dave Meko

% Whitens time series that are cols of array Y.
% Fits AR models (up to order 6) to the time series.
% Selects best model by AIC.
% Note that residual array E has the means of cols of Y added
%	back in, and that E has one or more leading zeros depending on the 
%	order model fit.  The user probably will want to lop off those 
%	leading years in E within the calling program before using E.

% Input arguments

%  Y (m1 x n1) the input time series
%  nhi - the maximum AR order to consider

% Output arguments

%   E (m1 x n1) the AR residuals, mean of Y added back in
%   K1 (1 x n1) order of fitted AR model
%   W (13 x n1) information on fit:

%	Cols 1:k1(i) - the AR coefficients
%	Cols 7:k1(i)+6 - the corresponding standard errors
%	Col 13       - ratios of residual to original variance

[m1,n1] = size(Y);

W=zeros(13,n1);
E=zeros(m1,n1);

mn=mean (Y);
Z=dtrend(Y);
NN=(1:nhi)';

for  i=1:n1
	z=Z(:,i);
	V=arxstruc(z,z,NN);
	nn=selstruc(V,'AIC');
	th=arx(z,nn);
	E(:,i)=resid(z,th)+ones(m1,1) * mn(i);
	k1=nn(1);
	W(1:k1,i)=th(3,1:k1)';
	W(7:k1+6,i)=sqrt(th(4,1:k1))';
	W(13,i)=th(1,1)/ (std(z) ^2);
	K1(i)=k1;
end


