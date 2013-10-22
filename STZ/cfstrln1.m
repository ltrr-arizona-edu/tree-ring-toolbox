function [cvx,a,b] = cfstrln1(yrv,yrvn,xvn,hwarn)
% cfstrln1:  fit straight line, either positive or negative slope
% 
% yrv -- years at which smoothed values are needed
% yrvn -- years at which original time series has data;  some might
%	be NaN
% xvn -- original time series;  some might be NaN
%
% cvx -- values of smooth curve at years yrv
%
% hwarn 1 if want warning window, 0 if not
%
% Fits data  to the equation
% g(t) = a + b*t
%   where t is the shifted time variable t = yrvn-yrvn(1)+1
%   In other words, t is same length as yrvn after 
%   dropping NaNs
%
%
% Any NaNs in yrvn and xvn should be internal, not at the ends of the
% series

tall = yrv-yrvn(1)+1; % Shifted time variable, all years, including
 % those with NaN data

% Remove NaNs
yrvn(isnan(xvn))=[];  % delete rows with xvn as NaN
xvn(isnan(xvn))=[];
nyr = length(xvn); % number of years for calibrating the fit

   
% Scale time series to begin at time t=1
t = yrvn-yrvn(1)+1;

% Fill data matrix, ones in col 1, shifted-year in col 2
W = [ones(nyr,1) t];

% Estimate parameters
c=regress(xvn,W);
%   Will return col vector, with constant in row 1, regression coef in row 2


% Estimated parameters
a=c(1);
b=c(2);

% Get smoothed curve at all time values, including those that might have
% been marked off with NaNs
cvx = a + b * tall;  


