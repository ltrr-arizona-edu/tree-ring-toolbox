function Y=armawhit(X,kopt,miss)
%Y=armawhit(X,kopt,miss)
%
% D Meko 9-2-95
%
% ARMA(1,1) residuals of time series from a time-series matrix
%
%
%***********************  IN ARGS ***************************
%
% X (mX x nX)r time-series matrix, year in col 1, variables in others
%		Missing values optionally filled with NaN or a specified code
%		(see kopt)
% kopt (1 x 1)i option for missing values
%	kopt(1)==1: NaN
%	kopt(1)==2: miss (see below)
% miss (1 x 1)r   missing value code
%
%********************** OUT ARGS ********************************
%
% y (mY x nY)r time series matrix or ARMA residuals
%
%
%********************** USER-WRITTEN FUNCTIONS CALLED -- NONE
%
%****************** NOTES ****************************************
%
% Time series in X may have variable time coverage, but no missing
% values internally (the series must be continuous)
%
% Either NaN or miss is read as a missing value, depending on kopt(1)
%
% Series in Y begin 1 year later than series in X because of loss of
% data in startup of modeling.  Y and X are same size, but leading 
% values of Y time series are set to missing value code


% Size and allocate
a=NaN;
[mX,nX]=size(X);
Y = a(ones(size(X)));


% Check X
if nX <2, error('X col-size must be at least 2'), end
% check that col 1 is consistent with being a "year"
yr = X(:,1);
dyr = diff(yr);
if ~all(dyr)
	error('X(:,1) is not a valid years column')
end


% Store data part of X
Z = X(:,2:nX);

% Check kopt(1), and make logical pointer to non-missing values
if kopt(1)==1
	b=NaN;
	Lvalid = ~isnan(Z);
elseif kopt==2
	b=miss
	Lvalid = Z~=miss;
else
	error('kopt(1) must be 1 or 2')
end


% Check that time series in Z are continuous
d = diff(Lvalid);
s1=sum(d== -1);
s2 = sum(d== 1);
if s1>1 | s2 >1, 
	error('One or more series in Z is not continuous')
end


% ARMA modeling
nsers = nX-1;
for n = 1:nsers;
	L1 = Lvalid(:,n); % pointer to non-missing data in Z
	z = Z(L1,n);
	u = z-mean(z); % convert time series to zero-mean
	th = armax(u,[1,1]);
	[e,r]=resid(u,th); % model residuals
	e = e+mean(z);
	plot(yr(L1),z,yr(L1),e);
	pause
	e(1)=NaN; % loss of startup value
	Y(L1,n+1)=e; 
end


Y(:,1)=yr;

