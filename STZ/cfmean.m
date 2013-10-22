function cvx = cfmean(yrv,yrvn,xvn)
% 
% Fit horizontal line through mean
%
% yrv -- years at which smoothed values are needed
% yrvn -- years at which original time series has data;  some might
%	be NaN
% xvn -- original time series;  some might be NaN

%
%
% cvx -- values of smooth curve at years yrv


  % Curve fitting
   yrvn(isnan(xvn))=[];  % delete rows with xvn as NaN
   xvn(isnan(xvn))=[];
   cvx=ones(length(yrv),1)*mean(xvn);
   

% End of file
