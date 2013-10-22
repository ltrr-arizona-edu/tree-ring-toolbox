function [Y,S,J]=gausmtx(X,p);
% gausmtx(X,p); gaussian filter the columns of a time series matrix
% [Y,S,J]=gausmtx(X,p);
% Last revised 1-26-02 
% 
% Filter the time series in a time series matrix by a gaussian filter with amplitude of frequency response 
% equal to 0.5 at a period of p years.  Requires the time series to be at least of length 2*p. 
%
%*** INPUT
%
% X (mX x nX)r  time series matrix to be filtered
% p (1 x 1)r  period at which amplitude of frequency repsons is 0.5
%
%*** OUTPUT 
%
% Y (mY x nY)r  filtered time series
% S (2 x nY)r   means (row 1) and standard devs (row 2) of the series in Y
% J (1 x ?)i   index to original rows of X indicating which series are the filtered series
%
%*** NOTES
%

%---- CHECK INPUT 

[mX,nX]=size(X);
if mX<20; 
    error('Row size of X must be at least 20');
end;

for n = 1:nX;
    x = X(:,n);
    % Check for internal NaN
    k=intnan(x);
    if k==1;
        error(['Internal NaN in col ' int2str(n) ' of X']);
    end;
    % Check for number of values in series
    ngood = sum(~isnan(x));
    if ngood<20;
        error(['Fewer than 20 nonNaN values in col ' int2str(n) ' of X']);
    end;
end;


% Computer filter length
b=wtsgaus(p);

t=(1:mX)';

Y=repmat(NaN,mX,nX);
Lpick=repmat(NaN,1,nX);



% FILTER SERIES

for n = 1:nX;
    x=X(:,n);
    L=~isnan(x);
    x=x(L);
    tx=t(L);
    mx = length(x);
    if mx<p;
        Lpick(n)=0;
    else;
        Lpick(n)=1;
        [y,ty]=filter1(x,tx,b,1);
        ny=length(y);
        irow = ty;
        Y(irow,n)=y;
    end;
end;

% UPDATE POINTER
Lpick=logical(Lpick);
J=find(Lpick);
if isempty(J);
    error('No series at least as long as filter period');
end;


% REDUCE MATRIX
Y=Y(:,J);

% STATISTICS
% Compute the means and statdard deviations of indices
S=[nanmean(Y);  nanstd(Y)];
