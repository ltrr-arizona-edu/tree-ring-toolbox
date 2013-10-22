function cvx = curvf(p,yrv,yrvn,xvn,nfit1)
% 
% p -- spline parameter
% yrv -- years at which smoothed values are needed
% yrvn -- years at which original time series has data;  some might
%	be NaN
% xvn -- original time series;  some might be NaN
% nfit1 -- curve-fit option
%  9 = spline
%  4 = mean
%
%
% cvx -- values of smooth curve at years yrv


  % Curve fitting
   yrvn(isnan(xvn))=[];  % delete rows with xvn as NaN
   xvn(isnan(xvn))=[];

  if nfit1==9;  % spline
   cvx=csaps(yrvn,xvn,p,yrv);
   cvx=cvx';
  elseif nfit1==4; % Horiz line through mean
   xvn(isnan(xvn))=[];
   cvx=ones(length(yrv),1)*mean(xvn);
  end


% End of file
