function [atilde,t] = pcareg1(L,W,y,G,tcrit)
%
% Principal Component Regression
%
% D Meko 8-23-95
%
%*****************8  IN ARGS *************************
%
% L (np x 1)r eigenvalues from PCA on input data
% W (mW x nW)r PCA amplitudes from PCA on input predictor data
% y (n x 1)r  predictand z-scores (might be PCA amps)
% G eigenvectors from PCA on X
% tcrit (1 x 1)r significance level for t-test on coefficients
%    (usually set tcrit to 0.95)
%
%**********************  OUT ARGS ******************
%
% atilde (p x 1)r   regression coefs for y on W

%*********************************  NOTES ***************
%
% For multi-predictand models, this function is called in a 
% loop to deal with each predictand separately.  Predictands
% are combined and predicted values transformed if necessary to
% original predictand variables in the calling program. 
%
% It is assumed that the predictors W are orthogonal to one
% another.  Thus the typical application would be for
% cols of W to be PC scores.  The predictand vector y might
% or might not be a PC-score series.
%
% Includes predictors only if estimator sig at 95% level
% Special case of none signif -- uses one with highest t 
%
% Method: First estimate model including all columns of W as 
% predictors.  Then delete predictors if their regression
% coefficients are too small as determined by a t-test.

% Full vector of estimators (regression coefficients)
ahat = inv(W'*W)*W'*y;

% Residual from full model
ehat = y - W*ahat;

n = length(y);
[mW,nW]=size(W);
p=nW;

t = (ahat .* sqrt(n*L)) / sqrt(ehat' * ehat/(n-p-1));
df = n-p-1; % degrees of freedom

tcrit = tinv(tcrit,df);
atilde= ahat;
[amax,imax]=max(abs(atilde));
L1 = abs(t)<tcrit;
s1 = sum(L1);
if s1>0 
	atilde(L1) = zeros(s1,1);
	if s1==p;
		atilde(imax)=amax;
	end
end


