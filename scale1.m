function yt=scale1(x,y,xt)
%
% Given time series data for x and y in a common period, and for
% x only in another period, estimate the data for y in the other period
% such that the standardized departures of y and x for the other period
% are identical.  
%
% Meko 2-20-97
%
%
%********************* IN *****************************************
%
% x (xm x 1)r source series for common period, to be used as predictor
% y (ym x 1)r corresponding segment of predictand series for common period
% xt (xtm x 1)r source series for target period
%
%*********************** OUT ***************************************
%
% yt (m1 x 1)r estimated y data for target period
%
%****************** METHOD *******************************************
%
% Compute the means and std devs of  x, y: xmean,ymean, xsd,ysd
% Compute the mean and std dev of xt:
%			xtmean,xtstd
%
% Compute standardized anomalies of x in target period using means from
% common period:
%
%    xsa= (xt-xmean)/xsd;
%
% Assume that y in the target period will have same stdzd anomalies from
% its common period mean:
%
%		yt = ymean + xsa*ysd;
%
%
% Algorithm assumes no NaN in x, y or xt
%
% Algorithm guarantees that in target period,
% standardized departures of yt are equivalent to standarized departures
% of xt, where the means and std devs for computing standardized departures
% are from the data in x and y
%
%*****************************************************************************

if any(isnan(x)) | any(isnan(y)) | any(isnan(xt));
	error('NaN found in x, y or xt');
end

% Compute means, standard devs
xmean=mean(x);
ymean=mean(y);
xtmean =mean(xt);
xsd=std(x);
ysd=std(y);
xtsd=std(xt);


% Compute standardized anomalies
xsa = (xt-xmean)/xsd;

% Transform to units scaled by mean, sdev of y
yt = ymean + xsa*ysd;

% Diagnostic debugging section -- not needed now
%figure(1);
%lent=length(xt);
%t=(1:lent)';
%plot(t,(xt-xmean)/xsd,t,(yt-ymean)/ysd);

