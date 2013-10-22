function [xbar,n]=meannan(X)
% [xbar,n]=meannan(X)
%
% Column mean of a matrix that might have NaNs as some elements
% NaNs ignored in computations
%
% Meko 1-17-97  Last revised 1-2-97
%
%*************** IN ARGS *******************************
%
% X (mX x nX)r tsm with mX observations on nX variables
%
%
%*************** OUT ARGS *******************************
%
% xbar (1 x nX)r means of the columns of X, ignoring NaNs
% n (1 x nX)i  sample size (number of observations) on which
%   mean is based
%
%*****************  NOTES ***********************************
%
% Trick to vectorizing is replacement of all NaNs with zeros while
% keeping track of the number of NaNs for each variable.  The zeros
% do not contribute to the column sum.
%
% To avoid "dividing by zero" errors, the sample size is changed to
% NaN from zero when when all elements of a column are NaN
%

a=NaN;

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

% Change column sample size to NaN if all elements of the column are NaN
L2=n==0; 
if any (L2);
	sum2=sum(L2);
	n(L2)=a(:,ones(sum2,1));
end


% Compute the means;  the zero values will not contribute to sums
sumx=sum(X);
xbar=sumx ./ n;




