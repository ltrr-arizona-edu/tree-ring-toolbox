function [gam,p90]=skew(x);

% Caution: invalid if sample size < 150

% Skewness coefficient and 90% CL, from Salas et al. (1980, p. 92)
% If series x is from a normal distribution, skewness is asymptotic
% normal with mean zero and variance 6/m, where m is the sample size.
% Equation for gam slightly different than eqn in WMO (1988, p 3 
% Appendix C) due to use of 

%  1.  1/m instead of m / (m-1)(m-2)   and
%  2.  1/m instead of 1/(m-1) in variance computation in denominator

%******   INPUT ARGUMENTS   *************************************
%
% x (m x n) time series
%
%********   OUTPUT ARGUMENTS   *********
%
% gam (m x 1)  skewness coeff for each col (variable) of x
% p90 (1 x 1)  90% of time, sample skew will fall within +-  p90
%	 of zero if time series is normally dist with skewness zero.
%
%******************************************************************



[m,n]=size(x);
xbar=mean(x);

if m==1;   % special case for row vector
	x=x';
	m=n;
	n=1;
end

xbar=mean(x);
dev=  x- xbar(ones(m,1),:);

sq= dev .^2;

num=(1/m) * sum(sq .* dev);

den1 = (1/m) * sum(sq);
den = sqrt(den1 .^3);

gam = num ./ den;

p90 = 1.645 * sqrt( 6 / m);




