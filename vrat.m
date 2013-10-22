function [F,dfb,dfw]=vrat(x,y)

% function to do analysis of variance test for significance of 
% apparent break in double mass curve

% From Bruce and Clark, 1969.  Introduction to Hydrometeorology. 
% Pergamon Press.  p. 162-165.

% Tests ratio of between sample variance to within sample variance.
% x and y are assumed to be col vectors holding data for two segments.
% F is the ratio of between to within sample variance, which should be
% tested with degrees of freedom 1 and N-2, where N is the total of
% the sample sizes of x and y.

% dummy trial data from bruce and clark
% x=[34.36 42.50 25.87 40.60 40.61 40.06 33.02 34.19 36.75 37.58]';
% y=[36.87 29.91 39.23 28.02 31.97 31.51 33.38 32.00 30.51 27.39]';

% [x y]

xsum= sum(x);
ysum=sum(y);

xsq=x' * x;  % sum of squares
ysq=y' * y;

T = xsum + ysum;  % total sum of all variables
N = length(x) + length(y);  % combined length of samples

CF= (T * T) ./ N;   % correction factor

TSQ = (xsq+ysq-CF);      % total sum of squares

dft= N-1;  % total deg of freedom

% compute between sample sos
BSSQ =((xsum * xsum) ./ length(x))+((ysum * ysum) ./ length (y))-CF;
dfb=1;

% compute within sample sos
WSSQ=TSQ - BSSQ;
dfw=dft-dfb;

% compute variances and variance ratios

BSV = BSSQ ./  dfb;

WSV = WSSQ ./ dfw;

F = BSV ./ WSV;



