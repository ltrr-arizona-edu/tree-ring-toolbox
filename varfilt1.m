function [wmean,wvar,yrw]=varfilt1(x,yr,m,kopt)
% [y,yry]=varfilt1(x,yr,m,kopt)
%
% Compute moving variance of sample mean
% D. Meko 4-24-98
%
%*******  INPUT ARGS
%
% x time series, col vect
% yr corresp years
% m   number of weights
% kopt options
%   kopt(1) which year of the m-year period to assign smoothed value
%		==1 central year  ==2  ending year
%
%**********  OUTPUT ARGS
%
% y   (?x1) moving variance (see Wilks, p. 123)
% yry (? x1) years for y

%**** Input check

% x and yr must be col vectors of same length
[mx,nx]=size(x);
[myr,nyr]=size(yr);
f=[nx~=1 nyr~=1 mx~=myr];
if any(f);
	clc
	error('x and yr must be col vectors of same length')
end

% yr must increment by 1 each "year"
if ~all(diff(yr)==1);
	clc
	error('yr does not have increment of 1')
end

% m must be positive integer, shorter than x
f=[m<=0  m>=length(x)  rem(m,1)~=0];
if any(f);
	clc
	error('m not postive integer, shorter than x')
end


% kopt must be row vector, length 1
[mkopt,nkopt]=size(kopt);
f=[mkopt~=1  nkopt~=1];
if any(f);
	clc
	error('kopt must be rv of length 1')
end;

% Compute how many m-year periods available in the time series
nper  = mx-m + 1;  

% Make a row vector of starting row increments; note that first period will 
% have increment of 0 from time series first year, second will have increment of
% 1, and so on to increment of nper-1;
inc1 = 0:(nper-1);
A = repmat(inc1,m,1); % dupe to a matrix

% Similar make an index of row indices to which the increment will be added
b = (1:m)';
B = repmat(b,1,nper);

% Make matrix of row indices to overlapping periods
C = A + B;  % size m x nper

% Pull data for periods, each in a col
W = x(C);

% Compute sample means and variances;
wmean = (mean(W))';
wvar = (var(W))';


i=(0:nper-1)';  % time vector for adding to years

% Compute "plotting" year vector
if kopt(1)==1; % centered moving average
	yr1=yr(1)+(m-1)/2;
elseif kopt(1)==2; % will want to plot smoothed series at end yr
	yr1=yr(1)+m-1;
end
yrw=yr1+i;

