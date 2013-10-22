function [xstd,n]=stdnan(X)
% [xbar,n]=stdnan(X)
%
% Column standard devs for matrix that might have NaNs as some elements
% NaNs ignored in computations
%
% Meko 1-17-97
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
%*********************  NOTES ******************************
%
% The input matrix X might contain some NaN elements.  
% The standard deviation is computed from the sum of squares
% of deviations from the mean. Substituting the column means
% for the NaNs leads to no contribution to the sums of
% squared deviations. 
% 
% Special case if all or all except 1 element of a column are 
% NaN.  Then the standard deviation cannot be computed, as 
% the denominator (n-1) would be 1 or 0 in computing the average
% squared deviation.  In this case, NaN is substituted for the
% sample size of 0 or 1.
%********************************************************

a=NaN;

% Make a copy of X
W=X;


% Size
[mX,nX]=size(X);
if mX==1 & nX>1; % X a row vector
	X=X';
	[mX,nX]=size(X);
	W=X;
end

% Compute sample size for each variable
L1=isnan(X);
sum1=sum(L1);
n= mX-sum1; % number of valid (not NaN) observations in each col of X

% Substitue zeros for any NaNs in X
if sum(sum(L1))~=0;
	X(L1)=zeros(sum(sum1),1);
end

% Change all sample sizes of zero to sample size NaN for computation of mean
L2=n==0;
if any(L2);
	sum2=sum(L2);
	n(L2)=a(:,ones(sum2,1));
end

% Compute the means;  the zero values will not contribute to sums
sumx=sum(X);
xbar=sumx ./ n;


% Substitute the col means for the NaNs in the copy of the
% original data 
XBAR=xbar(ones(mX,1),:);
W(L1)=XBAR(L1);

% Compute squared deviations from col means
D1=(W-XBAR) .^2;

% Sum squared deviations
sumd2 = sum(D1);

% Compute sample size for unbiased estimate of st dev
n2=n-1;

% Special case is n==0 or n==1. Then n2==-1 and n2==0, and 
% std dev cannot be computed.  Division by zero, in particular would
% give automatic warning message.  To avoid, want to generate NaN std dev
% if sample size is less than 2.  Already have changed zero elements
% of n to NaN.  Now need to do same for elements of n==1 (or elements of
% n2==0)
L3=n2==0;
sum3=sum(L3);
if sum3>0;
	n2(L3)=a(:,ones(sum3,1));
end

% Divide sum of squared deviations by valid sample size minus 1
xvar=sumd2 ./ n2;

% Square root to get standard dev from variance
xstd=sqrt(xvar);
