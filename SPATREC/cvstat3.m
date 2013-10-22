function [r2,mae,rmse,re]=cvstat3(y,yh,ymean)
% cvstat3: cross-validation stats for a single reconstructed time series
% [r2,mae,rmse,re]=cvstat3(y,yh,ymean)
%
% Meko 12-17-98
%
%********************* IN ARGS ************************
%
% y (my x 1)r actual predictand data, my years
% yh (same size as y)r cross-validation predictions for same
%		variables and years as covered by y
% ymean (my x 1)r calibration-period means for models used
%		to get the cross-validation predictions in yh
%
%
%********************** OUT ARGS *****************************
%
% r2 (1 x 1)r squared correlation coefficients between the 
%	cross-validation predictions and and actual data
% mae (1 x 1)r  mean absolute error of cross-validation predictions
% rmse (1 x 1)r root mean square error of cross-validation predictions 
% RE (1 x 1)r reduction-of-error statistic.  Based on mean-
%	-square-error of predictions and mse of hypothetical predictions
%	of calib-period mean for each year.

% Size and allocate
[my,ny]=size(y);
[myh,nyh]=size(yh);
if my~=myh | ny~=1 |  nyh~=1;
	error('y and yh must be column vectors of same length')
end
if size(ymean,1)~=myh | size(ymean,2)~=1;
   error('ymean must be col vector, same length as yh');
end

% Error quantities
e = y -yh; % residuals -- observed minus predicted
dm = y - ymean; % residuals -- hypothetical errors if the predictions had been
%    the calibration-period mean. These are also departures from mean
mae = mean(abs(e));  % mean absolute error of cross-val predictions
mse = mean (e .* e); % mean square error of cross-val predictions
rmse = sqrt(mse);
msemean = mean (dm .* dm); % mean square errors for prediction if calib
%   period mean had been used as the prediction in each year 
re = 1 - (mse ./ msemean);  % reduction of error statistic

% Squared correlation coefficient between actual and predicted for
% validation data
rr = corrcoef([yh y]);
r2=(rr(1,2)) .^2;
