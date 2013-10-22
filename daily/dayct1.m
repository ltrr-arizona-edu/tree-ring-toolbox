function [Y1,Y2]= dayct1(X,YRS,endday,c)
% dayct1: number of rainy days and total ppt in specified day window 
% CALL: [Y1,Y2]= dayct1(X,YRS,endday,c);
%
%**********************  IN ARGS  ***************************
%
% X (mX x nX) daily ppt matrix; no missing values, nX-2 stations
%		col 1: year
%`		col 2: day (1 to 366)
%		col 3:nX    ppt (inches)
% YRS (2 x 2)   first, last years of 
%		row 1: matrix X
%		row 2: desired summary period
% endday (1 x 2)   first, last days of summary window
%		Example:  [32 366] is Feb 1 through Dec 31
%			  [241 31] is previous year's Sept 1 (about) -Jan 31
% c (1 x 1)   threshold amount defining a "rainy" day;  typically c=0
%		Daily ppt must be greater than c for the day to be counted 'rainy'
%
%
%******************  OUT ARGS ***********************************
%
% Y1 (mY1 x nY1)   number of rainy days in day window (enddays)
%		at each station in each year of the specified analysis period 
%		col 1:  year
%		col 2-nY1:  values for each of nY1-1 stations
% 		rows:  each row is a year
% Y2 (mY2 x nY2)   like Y1, except total rainfall, in inches
%
%
%***********************  NOTES  ***********************************
%
% Previously name rainy.m when developed for Connie Woodhouse rainy day analysis
%
% Daily ppt matrix X should contain 366 days for each year, with possibly
% fewer days in the first and last year.  It is assumed that all missing 
% data have been filled in beforehand.  Day 60 (Feb 29) on non-leap years
% should have a missing-value code.  It makes no difference to this function
% what the code is, because the code is replaced with zeros for the count
% of number of rainy days and the total rainfall over the day window.
%
% Be careful to consider the time coverage by valid daily data in choosing
% the analysis period YRS(2,:).  For example, if valid daily data start in
% Jan 5, 1900 and end in Nov 28, 1993, the combination 
%  endday = [1 366], YRS = [1900 1993; 1900 1993] 
% gives bogus values for years 1900 and 1993.  If you want to have a day
% window 1-366 for this example, you can at best get an annual series
% for 1901-1993 (i.e., YRS=[1900 1993;  1901 1992]).
%
%


% Flesh out beginning and end of daily ppt matrix if first row is not 
% Jan 1 and last row is not Dec 31;  dummy values to be added are -99.
[mX,nX]=size(X);
daygo=X(1,2);
yrgo = X(1,1);
daysp=X(mX,2);
yrsp = X(mX,1);
mss = -99;
if daygo~=1 ; % first row of X is not for Jan 1
	daycv = (1:(daygo-1))';  % cv of days to splice at start of X
	ndd = length(daycv);
	yrcv = yrgo(ones(ndd,1),:);  % cv of duped first year
	MSgo = mss(ones(ndd,1),ones(nX-2,1));
	Xgo = [yrcv daycv MSgo];
	X=[Xgo; X];
end
if daysp~=366 ; % last row of X is not for Dec 31
	daycv = ((daysp+1):366)';  % cv of days to splice at end of X
	ndd = length(daycv);
	yrcv = yrsp(ones(ndd,1),:);  % cv of duped last year
	MSsp = mss(ones(ndd,1),ones(nX-2,1));
	Xsp = [yrcv daycv MSsp];
	X=[X; Xsp];
end

% Check data consistency
% Compatible row size of X and YRS?
[mX,nX]=size(X);
n1 = YRS(1,2) - YRS(1,1)+1; % specified # years in X
nr1 = 366 * n1;  % corresponding number of rows that should be in X
if mX ~= nr1
	error('Row size of X inconsistent with YRS(1,:)')
end

% Possible first year of day-windowed output, and first needed year
% of input data depend on whether day window crosses the calendar year
% boundary.
k1 = 0;  % indicator for day window crossing year boundary
if  endday(1) < endday(2);  % day window does not cross calendar year boundary
	first = YRS(2,1);  % will need data from this year to form first year's grouping
	goposs = YRS(1,1);  % first possible year for day-window data
else;  % crosses year boundary
	first = YRS(2,1) - 1;  % will need preceding years data
	goposs = YRS(1,1) + 1;  % earliest possible year for day-windowed series
	k1 =1 ;   % flag indicating that day window crosses year boundary
end
	
% Specified years for analysis period compatible with year coverage of X?
if YRS(2,1) < goposs
	error('YRS(2,1) start year for analyis imposssible for this enday and X')
end
if  YRS(2,2) > YRS(1,2)
	error ('YRS(2,2) later than coverage by X')
end


nyrs = YRS(2,2)-YRS(2,1)+1;  % # of years of output series
yr = (YRS(2,1):YRS(2,2))';  % cv of output years
nstns = nX - 2;  % number of stations

%count number of missing value (9999) in the daily sounding data
L=X(:,3)==9999;
sum(L);

%substitute missing value with NaN
X(L,3)=NaN;

%compute long term daily mean and standard deviation

Xmn=repmat(NaN,366,2);
Xstd=repmat(NaN,366,2);

Xmn(:,1)=(1:366)';
Xstd(:,1)=(1:366)';

for i=1:366
   L5=X(:,2)==i;
   XX1=X(L5,3);
   xmn=nanmean(XX1);
   Xmn(i,2)=xmn;
   xstd=nanstd(XX1);
   Xstd(i,2)=xstd;
end



% Compute number of days in day window, and form pointer to days
if k1==0;  % not cross year boundary
   L2 = X(:,2) >= endday(1) & X(:,2)<= endday(2);
   L6 = Xmn(:,1) >= endday(1) & Xmn(:,1) <= endday(2);
	ndays = endday(2) - endday(1) +1;  % # days in day window
else
   L2=(X(:,2)>=endday(1) & X(:,2)<=366)|(X(:,2)>=1 & X(:,2)<= endday(2));
   L6=(Xmn(:,1)>=endday(1) & Xmn(:,1)<=366)|(Xmn(:,1)>=1 & Xmn(:,1)<= endday(2));
	ndays = (366 - endday(1)+1) + endday(2);
end


% Make pointer to years of X needed for analysis period
L1 = X(:,1) >= first & X(:,1) <= YRS(2,2);


% Initialize output matrices and put year in col 1
na  = NaN;
Y1 = na(ones(nyrs,1),ones(nstns+1,1));
Y2 = na(ones(nyrs,1),ones(nstns+1,1));
Y1(:,1) = yr;
Y2(:,1) = yr;


%  Get needed rows of X, Xmn and Xstd
X1 = X(L1 & L2,:);
Xmn1=Xmn(L6,2);
Xstd1=Xstd(L6,2);

% Truncate leading and trailingif window crosses year boundary;
if k1==1;
	X1(1:endday(2),:)= [];
	[mX1,nX1]=size(X1);
	nout = 366 - endday(1) + 1;  % # trailing days to drop off
	nn = mX1 - nout;
	X1=X1(1:nn,:);
end


% X1 should now have only the specified days in the day window, for the
% years needed to produce the output for the period YRS(2,:)
[mX1,nX1]=size(X1);
nr2 = ndays * (YRS(2,2) - YRS(2,1) + 1);  % expected number of rows in X1
if nr2 ~= mX1
	error('Row size of X1 incorrect')
end


% LEAP YEAR ADJUSTMENT:  If day 60 (Feb 29) is in the day window, replace
% the existing daily ppt value for non-leap years (a missing-value code)
% with 0.  Then we can work with constant-size matrices for each year,
% making vectorized computations possible.  The substitution does not 
% affect the total rainfall or the the number of rainy days.  If you
% want number of non-rainy days using ndays-# of rainy days,  however, 
% you will need to increase value by 1 for leap years
L3 = (rem(X1(:,1),4) ~= 0) & X1(:,2)==60;
zz = zeros(1,nstns);
if any(L3)
	zz = zz(ones(sum(L3),1),:);
	X1 (L3,3:nstns+2) =  zz;
end

% Check for missing value codes in your culled out daily data.
% If any such codes, you must not have considered valid data
% coverage properly when specifying the analysis years and
% day window
LM = X1==9999;
if any(any(LM))
  error('9999 found in selected daily segment')
end

% Compute number of low (trough) days and number of high(ridge) days in day window
for n = 1:nstns ;   % loop over stations
	j = n+2;  % pointer to column in X1
	x2 = X1(:,j);  %  get a stations' daily data
   x2 = reshape(x2,ndays,nyrs);  % make into a matrix convenient for summing
   c1=Xmn1+Xstd1; % one standard deviation above long term daily mean -- for ridge 
   c2=Xmn1-Xstd1; % one standard deviation below long term daily mean --for trough
   % Dupe col vectors c1 and c2 to matrices same size as x2
   C1=repmat(c1,1,nyrs);
   C2=repmat(c2,1,nyrs);
	L4 = x2 > C1;   % indicator -- 1 if presure heights above the critical point
	y1 = (sum(L4))';  % number of high(ridge) days
	L5 = x2 < C2;   % indicator -- 1 if pressure heights below the critical point
	y2 = (sum(L5))';  % number of low (trough) days
	
	Y1(:,n+1) = y1 ;  % put in slot in output matrix
	Y2(:,n+1) = y2;
end


