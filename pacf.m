function [phi,SE2] = pacf(z,nlags)
% pacf: partial autocorrelation function and its standard error
% [phi,SE2] = pacf(z,nlags);
% Last revised 9-2-99
%
% Compute partial acf and its standard error for a single time series
%
%*** INPUT ARGS
%
% z  (? x 1) time series
% nlags (1 x 1)   number of lags for pacf
%
%*** OUTPUT ARGS
%
% phi (1 x nlags) partial autocorrelation, lag 1 thru nlags
% SE2 (1 x nlags) two-standard errors for pacf, lags 1 thru nlags
%
%*** REFERENCES
%
% O.D. Anderson, 1976, Time series analysis and forecasting.  
%  Butterworths, London and Boston, p. 9-10
%
%*** UW FUNCTIONS CALLED -- none
%*** TOOLBOXES NEEDED 
% system identification
%
%*** NOTES


phi= ones(1,nlags);  % initialize length of return vector

z=dtrend(z);     % subtract mean from z
[m,n]=size(z);
P=zeros(m,m);  % preallocate an array to hold the autoc mtx

r=covf(z,m)';  % autocovariance function, lags 0 thru m
r1=(r ./r(1)') ;  % autocorrelation function

s=r1(m:-1:2,1);  % reversed array r1
g=[s' 1 r1(2:m)']; % row vector, first m-1 elements are autocorrs
	% for lags m-1 thru 1;  middle element is "1";  last m-1
	%elements are autocorrs for lags 1 thru m-1

% Build the autocorrelation matrix

for i= [1:1:m]
	P(i,:) = g(m-i+1:2*m-i);
end

% loop to compute partial autocorrelations

phis(1)=r1(2);  % first pac is equal to first autocorr coef

for l=[2:1:nlags]
	psub=P(1:l,1:l);
	pstar=[psub(:,1:l-1)  r1(2:l+1)];
	phis(l) = det(pstar)/det(psub);
end

phi(1:nlags)= phis(1:nlags);

se2 = 2.0 / sqrt(m);    %  Two std errors, as function of sample size
SE2=se2(:,ones(nlags,1)); % Convert to rv
