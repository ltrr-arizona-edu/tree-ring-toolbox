function dummy = acfpacf(xdata)
%   acfpacf.m ...  compute and plot acf and pacf for time series
%   in cols of a matrix

%  3-21-92 by D. Meko:  for demonstration screen plots of acf and
%      pacf, with 2-se bars

%******  USER-WRITTEN FUNCTIONS NEEDED
%
% acf.m   --- to get acf for a vector time series
% pacf.m  --- to get pacf for a vector time series


%********  INPUT ARG
%
%  xdata (? x ?)?  screen-prompted for name of data array, 
%        may or may not have year in col 1

%*******  OUTPUT 
%
% No data output, only plots of acf and pacf with 2-se bars.  See function
% definitions of acf.m and pacf.m for description of methods.

%******  BEGIN CODE

k1=[];
k1=input('IS THERE A YEAR COLUMN IN DATA MATRIX? Y/N [N]','s');
if isempty(k1), k1='N';  end;
if k1=='Y',
	xdata(:,1)=[];
end


[m,n]=size(xdata);  % Find number of rows and cols in xdata.
nlags=25;  % Call for values of acf and pacf at lags 1 to 25
L=0:nlags;  % to be used as x-axis on plots

for i=1:n;  % Loop for each time series
	clg
	subplot
	z=xdata(:,i);  % cv holding current series to be analyzed
	
	[r,sr,r95]=acf(z,nlags);  % Compute autocorr function, lags 1-25
	%     sr holds two-standard errors (see function definition)
	[phi,sp]=pacf(z,nlags);  % Compute partial acf, lags 1-25

	subplot(211)
	plot(L,[1 r],'*',1:nlags,sr,'--',1:nlags,-sr,'--')
	title(['SERIES # ',int2str(i),'  ACF AND 2 LARGE-LAG SEs']);
	pause
	plot(L,[1 phi],'*',1:nlags,sp,'--',1:nlags,-sp,'--')
	title(['SERIES # ',int2str(i),' PACF AND 2 SEs']);
	pause

end


