function [tau,point95]=mannkend(x,y)

% Mann-Kendall test for inhomogeneity (WMO, 1966, p. 65)
% Pgmd 2-6-92 by D. Meko
% Modified to correct bug 2-22-92

%***  INPUT ARGUMENTS
%
% x (cv)  a precipitation record
% y (cv)  another precipitation record, same years as for x


%******  OUTPUT ARGUMENTS

% tau (1 x 1) test statistic
% point95     95% probability point for dist of tau

%*******  NOTE:  The handling of zero seasonal ppt assumes data is
%		passed to this function in units of hundredths of inches.


N=length(x);
M=length(y);


if N~=M
	error('SERIES UNEQUAL IN LENGTH')
end

% Do not want infinite ratio in x / y.
% If both zero value for both x and y, replace each with a small
% non-negative number.  Resulting ratio will be 1.
% If zero value in denom only,  replace zero values in denom series
% with a small non-neg value.  If x, y are in hundredths of inches
% for example, replace zero values with 1, which equals 0.01 in.
lyy= y==0 & x==0;
if any(lyy);
	y(lyy)=1.0;
	x(lyy)=1.0;
end
ly = y==0  ;  % point to zero values in denominator series
if any(ly)
	y(ly)=1.0;  % Change this line depending on units of ppt in x, y
end

% Form ratio series
	x= x ./ y;

for i = 1:N-1
	a = x(i);
	n(i) = sum(a<x(i+1:N));
end

P=sum(n);
tau = ((4*P) / (N*(N-1))) -1;
point95 = 1.96 * sqrt((4*N+10)/(9*N*(N-1))); % 1.96 is 95% prob point
		% of standard normal dist for 2-tailed test
