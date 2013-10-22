function H=hinge(X)
% H=hinge(X)
% 
% Median, fourths, extremes, for columns of matrix X
% Computed following Hoaglin, Mosteller and Tukey (1983):
% Understanding Robust and Exploratory Data Analysis. John
% Wiley & Sons, Inc.
%
% D Meko 2-28-96
%
%**************** IN ARGS ******************
%
% X (mX x nX)r  data matrix.  Each column might be, say,
%	an overlapping 30-year period of a ring-width series
%
%******************** OUT ARGS ****************
%
% H (5 x nH)r    lower extreme, lower fourth, median, upper
%	fourth, upper extreme for each of the columns of X.  Note
%  that nH==nX
%
%**************** USE *********************************
%
% Need fourths for power transformations described by 
% Hoaglin et al. (1983).  Say have a ring-width series of 
% 350 years. Say decide to compute fourths and medians for 
% overlapping 30-yr segments.  Would first call a function
% like pullseg1.m 
% to segment the original time series so that X would be
% a matrix whose columns are the segments of the original 
% series.  Then would call hinge(X) for the fourths, etc.
% Then would call some other function to plot functions of
% the medians and fourths to arrive at a suitable power
% transformation.


% Check inputs, size and allocate
[mX,nX]=size(X);
a=NaN;
H = a(ones(5,1),ones(nX,1));
if mX==1, 
	error('X must be column vector or matrix, not a row vector')
end

n = mX; % batch size,for simple reference
% Sort cols of X from smallest to largest values
S=sort(X);


% Medians
r = rem(n,2); % zero if n even, 1 if n odd
if r==0; % if number of observations in batch is even
	k=n/2;
	dm = k+0.5; % depth of the median
	m=0.5 * (S(k,:) + S(k+1,:)); % rv of medians
else; % odd number of observations in batch
	dm=(n+1)/2; % depth of the median
	m = S(dm,:);
end
H(3,:)=m;


% Extremes
H(1,:) = S(1,:);
H(5,:)= S(n,:);


% Fourths
dL = (floor(dm)+1)/2;  % depth of fourths, also is rank of
		%		lower fourth
dU = n + 1 - dL; % rank of upper fourth
if rem(dL,2)~=0; % if depth of fourths non integer
	fL = (S(floor(dL),:) + S(ceil(dL),:))/2;    % lower fourth
	fU = (S(floor(dU),:) + S(ceil(dU),:))/2;   % upper fourth
else; % depth of fourths is integer
	fL = S(dL,:);
	fU = S(dU,:);
end
H(2,:)=fL;
H(4,:)=fU;
