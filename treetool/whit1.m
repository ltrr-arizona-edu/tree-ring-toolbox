function [e,k1,vrat,arcs] = whit1(y,nhi,k2)
% whit1:  fit AR model to a time series, returning model information and residuals
% [e,k1,vrat,arcs] = whit1(y,nhi,k2);
% Last revised 7-12-99
%
% You specify the highest order AR model to consider.  Models up to that order are
% fit, and the Akaike Information Criterion (AIC) is used to pick the best model.
%
%*** IN **********************
%
% y (my x 1)r  time series, vector; NaNs not allowed
% nhi (1 x 1)i  highest AR order to consider (if k2==1), or the only
%   AR order model to try fitting (if k2=2)
% k2 (1 x 1)i  option for order selections
%   ==1 fit models of order 1 to nhi
%   ==2 fit model of order nhi only
%
%*** OUT **************************
%
%  e (my x 1)r AR residuals, with mean added back in; or, if null model, the
%		original series y
%  k1 - the order of the AR model deemed best by the AIC, or 0 if null model.
%  vrat - the ratio of variance of AR residuals to variance of original time series y
%  arcs the ar coefs and their two-standard errors in a two-row array; row 1 has
%   the coefficients; row two has 2*standard error of the coefficients 
%
%*** REFERENCES
% 
% The loss function and AIC are discussed in Ljung, L. 1995. System identification 
% toolbox; for use with MATLAB, The MathWorks, Inc., p. 3-46.
% 
%*** UW FUNCTIONS CALLED -- none
%*** TOOLBOXES NEEDED -- system identification
%
%*** NOTES
%
% vrat: the variance ratio is computed over years k1+1 to n1, where n1 is length of y
%
% null model (k1==0): 
% The "Loss function" V is the normalized sos of prediction errors, or the 
% variance of the AR residuals.
% A "null model" loss function is defined as the mean sos of the
% variable z (original data with mean subtracted) over the years nhi+1 to the
% last year of z.  If the null-model loss function is no higher than the V
% for any estimated model, a "null model" is assumed.  This means no AR model is
% justified, and k1==1 is returned.

% 

mn=mean(y);
z=dtrend(y);
n1=length(y);

if k2==1;   %  mode 1, fit up to order nhi
	NN=(1:nhi)';  % orders to try
else;  % mode 2, fit order nhi only
	NN=nhi;
end

V=arxstruc(z,z,NN);
[nn,vmod]=selstruc(V,'AIC');

th=arx(z,nn);
[e,r]=resid(z,th);

k1=nn(1);


nsub=n1-nhi;
vmodnull = log(((nsub-1)/nsub) * (std(z(nhi+1:n1))) ^ 2);

vrat= (std(e(k1+1:n1)) ^2) / (std(z(k1+1:n1)) ^2);

e=e+mn;

thsub=diag(th(4:3+k1,1:k1));  % variances of est params
stderr= 2*sqrt(thsub);
arcs=[th(3,1:k1); stderr'];

% Check that null model (ar order =0) is not better than any of others
% considered.

L1=vmodnull(:,ones(length(vmod),1)) <=  vmod(1,:);
if all(L1);   %  Loss function for null model is smaller than for
		% any other model, so do not whiten.  Pass the orig
		% series back as the "whitened" series, and tell
		% that model is null

	k1=0;
	e=y;
	vrat=1.0;
	arcs=[0 0]';	
end


