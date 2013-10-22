function [R,pval]=portmant(r,N,p,q,K)
% portmant: portmanteau lack-of-fit test on ARMA model residuals
% [R,pval]=portmant(r,N,p,q,K);
% Last revised 2-15-98
%
%***  IN
%
% r (1 x mr)r  autocorrelation function at lags 1 through mr
% N (1 x 1)i  sample size (number of obs) of the time series r based on
% p (1 x 1)i   AR order of model
% q (1 x 1)i   MA order of model
% K (1 x 1)i  <20> number of lags on which the portmanteau statistic is to be computed
%
%*** OUT 
%
% R (1 x 1)i  computed portmanteau statistic
% pval (1 x 1)r  p-value for chi-squared test.  A small pval indicates
%   reject null hypothesis that first K autocorrelation coefficients as a 
%   group come from a series of shocks that is random (adequate model)
%   For example, if pval is 0.009, would reject (at the 0.01 level) 
%   the null hypothesis that the model is adequate
%
%*** REFERENCES 
% Anderson, 1976.  Time series analysis and forecasting. 
% Butterworths, Boston,  182 pp.  Page 84 describes the portmanteau statistic.
%
%*** UW FUNCTIONS CALLED -- none
%*** TOOLBOXES NEEDED -- statistics
%
%************* NOTES
%
% User previously has fit and ARMA(p,q) model to a time series of length N, and
% has generated the model residuals and computed their acf   r.
% A significant R indicates model inadequacy
%
% This is a low-powered test, in that if the statistic proves not significant, 
% not much support can be given to the model.

%
if nargin == 4;
   K = 20;
end

if length(r) <K;
   error('r must be at least length K');
end

   
% ********* COMPUTE STATISTIC

ruse = r(1:K); % use first K ac coeffs in computing R
R = N * sum(ruse .* ruse);   % eqn 9.10 in Anderson, 1976

%************ COMPUTE DEGREES OF FREEDOM and PROBABILITY LEVEL
df = K - p -q;
pval = 1-chi2cdf(R,df);



