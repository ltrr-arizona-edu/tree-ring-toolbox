function [s,n]=sumnan(X)
% [s,n]=sumnan(X)
%
% Column sum of a matrix that might have NaNs as some elements
% NaNs ignored in computations
%
% Meko 1-29-97
%
%*************** IN ARGS *******************************
%
% X (mX x nX)r tsm with mX observations on nX variables
%
%
%*************** OUT ARGS *******************************
%
% s (1 x nX)r sums of the columns of X, ignoring NaNs
% n (1 x nX)i  sample size (number of observations) on which
%   sum is based
%
%*****************************************************

% Size
[mX,nX]=size(X);
if mX==1 & nX>1; % X a row vector
	X=X';
	[mX,nX]=size(X);
end

% Compute sample size for each variable
L1=isnan(X);
sum1=sum(L1);
n= mX-sum1; % number of valid (not NaN) observations in each col of X

% Substitue zeros for the NaNs in X
X(L1)=zeros(sum(sum1),1);

% Compute the sums;  the zero values will not contribute to sums
s=sum(X);





