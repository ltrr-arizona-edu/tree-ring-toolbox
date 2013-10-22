function [r,SE2,r95]=acf(x,nlags)
% acf: autocorrelation function and approximate 95% confidence bands
% [r,SE2,r95]=acf(x,nlags);
% Last revised  2-2-98
%
% Computes autocorrelation function and related information for a time series.  
%
%*** INPUT ARGUMENTS
%
% x (m x 1) time series, length m
% nlags (1 x 1)i number of lags to compute acf to
%
%*** OUTPUT ARGUMENTS
%
% r (1 x nlags)r   acf at lags 1 to nlags
% SE2 (1 x nlags)r two times the large-lag standard error of r
% r95 (1 x 1)r value of first-order autocorrelation coef (r(1))
%		significant at 95% level in one-tailed test (WMO method)
%
%*** REFERENCES
%
% Large-lag standard error after Box, G.E.P., and Jenkins, G.M., 1976, Time series 
% analysis: forecasting and control: San Francisco, Holden Day.  
%
% Probability point for one-tailed test of first-order autocorrelation coefficient 
% after World Meterorological Organization, 1966, Technical Note No. 79: Climatic 
% Change, WMO-No, 195.  TP.100, Geneva, 80 pp.
%
%
%*** UW FUNCTIONS CALLED -- none
%*** TOOLBOXES NEEDED - none
%
%*** NOTES  ****************************************
%
% Large-lag standard error from p. 35, Box and Jenkins (1976)
% Prob point for one-tailed test from p. 60,  WMO (1966)
%
% Modifications:
%   Original function written in 1991
%   Mod 10/92 to give results for lags 1-M instead of 0-M

[m,n]=size(x);
if n~=1
	error('FUNCTION ACF REQUIRES VECTOR TIME SERIES')
end

c=covf((x - mean(x)),nlags+1);  % autocovariance function, lags 0 to nlags+1
r=c(2:nlags+1) / c(1);  % divide by variance

r95 = (-1 + 1.645*sqrt(m-2)) / (m-1);  % for one-tailed test, 95% signif.


%************   LARGE-LAG STANDARD ERROR   ******************


for i=1:nlags
	if i==1  ;   %  sum term zero for first-order r
		sum1 = 0;  
	else
		rsub=r(1:i-1); 
		rsubsq=rsub .* rsub;
		sum1=2*sum(rsubsq);
	end

	var=(1/m)  *  (1.0 + sum1);  % variance of ac coef
	SE2(i) = 2.0 * sqrt(var);  % two standard errors
end

	

