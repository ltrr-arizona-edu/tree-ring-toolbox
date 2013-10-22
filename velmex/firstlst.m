function T = firstlst(X,yrX)
% firstlst:  find first and last years of valid data in a time series matrix
% T = firstlst(X,yrX);
% Last revised 2005-8-22
%
% Find first and last years of valid data in a time series matrix
%
%*** IN
%
% X (mX x nX)r time series matrix
% yrX (mX x 1)i year (or time) vector for X
%
%
%*** OUT
%
% T (nX x 2)i  first and last years (or observation numbers) of series in X
%
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% The input matrix X has time series with possibly differing time coverage by the variables in columns. 
% T lists the first and last year of non-NaN data for each series

[mX,nX]=size(X);
if ~(isvector(yrX) && length(yrX)==mX);
    error('yrX must be col vector same row-size as X');
end;

d = diff(yrX);
if ~all(d>=1);
    error('yrX must be monotonic increasing');
end;

D = repmat(yrX,1,nX); % yrX duped to matrix same size as X
L = isnan(X); 
D(L)=NaN;
T = [nanmin(D); nanmax(D)]';





