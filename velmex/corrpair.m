function [r,n]=corrpair(X,Y)
% corrpair: Pearson correlation between pairs of columns in two matrices that may contain some NaNs
% R=corrpair(x,m);
% Last revised 2006-9-2
%
% Pearson correlation between pairs of columns in two matrices that may contain some NaNs. 
% Element r(1) is the correlation between col 1 of X and col 1 of Y;  element r(2) is the correlation 
% between col 2 of X and col 2 of Y; etc. For any pair of series, only those observation with valid data for
% both series are used for computing the quantities needed for correlations (e.g., standard deviations,
% mean cross-product, sample size)
%
%*** INPUT
%
% X (mX x nX)r  first matrix
% Y (mY x nY)r second matrix
%
%*** OUTPUT
%
% r(1 x nX)r  correlation coefficients
% n(1 x nX)r  sample size for correlation coefficients
%
%*** REFERENCES
%
% Wilks, D.S., 1995, Statistical methods in the atmospheric sciences: Academic Press (p. 46)
%
%*** UW FUNCTIONS CALLED -- NONE
%
%*** TOOLBOXES NEEDED
%
% statistics
%
%*** NOTES
%
% X and Y must be same size


%--- CHECK INPUT

[mX,nX]=size(X);
[mY,nY]=size(Y);
if mX~=mY | nX~=nY;
    error('X and Y must be same size');
end;


%--- TALLY SAMPLE SIZE FOR EACH PAIR (NOT INCLUDING NANS)

L = ~isnan(X) & ~isnan(Y);
n = sum(L); % rv of sample sizes
if any(n<3);
    error('Some pairs of series have fewer than three valid (non-NaN) observations for computing correlation');
end
N = repmat(n,mX,1); % duped rows to matrix of sample size


%--- ADJUST SERIES SO EACH IN PAIR HAS NAN WHEN ONE HAS NAN

X(~L)=NaN;
Y(~L)=NaN;


%--- COMPUTE MEANS OF ADJUSTED SERIES
meanx = nanmean(X) ; % rv of means of X
meany = nanmean(Y); % rv of means of Y 
MX = repmat(meanx,mX,1);
MY = repmat(meany,mY,1);


%--- COMPuTE DEPARTURES FROM (censored) MEANS

DX = X-MX;
DY = Y-MY;


%--- COMPUTE STANDARD DEVIATIONS

sx = nanstd(X);  % rv of standard devs for X
sy = nanstd(Y);;  % rv of standard devs for Y


%--- COMPUTE SAMPLE COVARIANCES

pd =   (nansum(DX .* DY)) ./  (n-1);


%--- COMPUTE CORRELATIONS

r = pd ./ (sx .* sy);
