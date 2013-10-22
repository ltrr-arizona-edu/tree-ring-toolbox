function Y=truncass(X,ptyr,cutyr)
% truncass:  truncate start years of chrons depending on sample size
% CALL: Y=truncass(X,ptyr);
%
%******************  IN **************************************
%
% X (mX x nX)r  tsm of standard tree-ring chrons, with year in col 1
% ptyr (m1 x 1)i  'point year' or 'ass year' representing first year
%   with acceptable sample size in series in X.  Note that m1==nX-1
% cutyr (1 x 1)i  cutoff year. Do not change any real data to NaN for this
%   year on no matter what the ptyr is
%
%******************** OUT **************
%
% Y (mY x nY)r tsm of truncated chrons. Y is same as X, except data for years
% before point years are changed to NaN
%
%
%*********** NOTES
%
% Written for cleaning up series before using arwhite2.m in \af1\sacrflow study

[mX,nX]=size(X);
[m1,n1] = size(ptyr);


Y=X;  % initialize Y; will hold the truncated versio nof X

if n1~=1;
   error('ptyr must be cv');
end

% Adjust ptyr for cutyear
L1 = ptyr>cutyr;
ptyr(L1)=cutyr;


nsers = nX -1; % number of time series in X
if m1~=nsers;
   error('row size of ptyr should be 1 less than col size of X');
end

% store year col of X
yr = X(:,1);

% Remove year col
X(:,1)=[];
Y(:,1)=[];

% Loop over series
for n = 1:nsers;
   x=X(:,n);
   L=yr<ptyr(n);
   x(L)=NaN;
   Y(:,n)=x;
   plot(yr,X(:,n),yr,Y(:,n));
   title(['Series ' int2str(n)]);
   grid;
   pause(1);
end

% Put year col back
Y = [yr Y];









