function [R2,RE]=cvstat1(Y,Yh,M)
% [R2,RE]=cvstat1(Y,Yh,M)
%
% On output from spatrec1.m, compute crossvalidation statistics for
% single-site regression models.
%
% D Meko 12-24-95
%
%********************* IN ARGS ************************
%
% Y (mY x nY)r actual predictand data, mY years, nY sites
% Yh (same size as Y)r crossvalidation predictions for same
%		variables and years as covered by Y
% M (same size as Y)r calibration-period means for models used
%		to get the predictions in Yh.
%
%
%********************** OUT ARGS *****************************
%
% R2 (nY x 1)r squared correlation coefficients between the 
%	predicted and actual data (Yh and Y).
% RE (nY x 1)r reduction-of-error statistic.  Based on mean-
%	-square-error of predictions and mse of hypothetical predictions
%	of calib-period mean for each year.

% Size and allocate
[mY,nY]=size(Y);
[mYh,nYh]=size(Yh);
if mY~=mYh | nY~=nYh,
	error('Y and Yh must be same size')
end
a=NaN;
R2 = a(ones(nY,1),:);
RE=a(ones(nY,1),:);

% Reduction of Error Statistic
E = Y -Yh; % residuals -- observed minus predicted
DM = Y - M; % residual -- observed minus calib-period mean
MSE = mean (E .* E); % rv of mean square errors for prediction
MSM = mean (DM .* DM); % rv of mean square errors for prediction of 
		% calib-period mean
RE = 1 - (MSE ./ MSM);  % reduction of error statistic
RE=RE'; % rv to cv

% Squared correlation coefficient between actual and predicted for
% validation data
for n = 1:nY;
	yh = Yh(:,n);
	y = Y (:,n);
	rr = corrcoef([yh y]);
	R2(n)=(rr(1,2)) .^2;
end
