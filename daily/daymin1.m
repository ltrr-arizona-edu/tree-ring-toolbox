function [yr,q,T]=daymin1(X,n,daywind,kopt)
% daymin1:  compute annual time series of n-day minimum flows in a seasonal time window
% [yr,q,T]=daymin1(X,n,daywind,kopt);
% Last revised 2-15-02
%
% Compute annual time series of n-day minimum flows in a seasonal time window.  Typical use is after 
% generating the tsm of daily flows (e.g., with dayusgs1.m).  A time window is defined for each year. 
% For a given year, running n-day means of daily flow are computed, and the minumum running meab is
% recorded, along with the end day of the lowest-flow period.  The process is repeated for every year.  
% The time window can cross the break in calendar years.  Missing data are acceptable, as long as there
% is at least one n-day period with data in the time window.  If not, the n-day-minimum for the year is marked
% as NaN.
%
%*** INPUT
%
% X (? x 5)r  daily flow time series matrix -- see daycon2.m or dayusgs1.m  for format
% n (1 x 1)i  length of n-day period
% daywind (1 x 2)i window defining period for which n-day minimum is 
%    selected.  daywind(1), daywind(2) are endpoints of window in 366-day
%    terms.  For example daywind==[1 366] means the n-day minimum whether
%    it occurs anytime in period Jan 1 to Dec 31.  Or [32 60] means the n-day
%    minimum occuring in February.  For water year, use [275 274].  This 
%    specifies Oct1 thru Sep30 (see notes)
% kopt(1 x 1)i  unused option
%
%
%*** OUTPUT
%
% yr (? x 1)i  year vector for the n-day=minimum flows.  The year in yr is the year containing the
%   ending day of the time window (see notes)
% q (? x 1)i  n-day-minimum flow for each year in vector yr (cfs)
% T (? x 4)i  reference time matrix  for the the data vector q.  The cols are:
%   1 year, identical to yr
%   2,3  month and day-of-month of the last day of the n-day period with the minimum n-day flow for the year
%   4 - the Julian day corresponding to the month and day-of-month in cols 2,3 (see notes)
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED
%
% ndayma1 -- to compute moving averages
%
%
%*** NOTES
%
% daywind.  If daywind(1)>daywind2, window is understood to cross the calendar
%    year boundary (Dec 31), and "year" associated with a value of q is the 
%    year of the endpoint of the window.  Thus, if daywind is [275 274], the
%    min flow for Oct 1 1950 to Sept 30 1951 is linked with "year" 1951. The 
%    n-day minimum must END in the "year" to be listed with that year. Thus
%    for the water year, if the 5-day minimum for WY 1951 happens to comprise
%    Sep 27 1950 - Oct 1 1951,   that is the 5-day minimum for 1951.  
%
% yr. If the window crosses the calendar year, yr is the ending year of the seasonal window.  For example
%   if the window is set as [330 30], roughly Dec through Jan, the year in yr is that containing January
% T.  Example of a row of T:
%     1905  2  12 43
%     means that for the 1905 "year", the n-day-minimum flow was in the period ending Feb 12, or Julian day 43


% Set months strings
nmsmon={'Jan','Feb','Mar','Apr','May','June','July','Aug','Sept','Oct','Nov','Dec'};


%******  Compute number of years in the search window & check consistency
if any(daywind>366) | any(daywind<1);
   error('daywind elements must be between 1 and 366');
end;
if daywind(2)>=daywind(1);
   nwind = daywind(2)-daywind(1)+1;
else;
   nwind = daywind(2) +   (366-daywind(1)+1);
end;
if n>nwind; % cannot compute n-day low if annual period shorter than n days
   error('n greater than width of day window');
end;


%***** Compute number of years with fully available windows and if necessary
% lop off first year so that all search windows complete

% How many occurences of ending day?
L1 = X(:,4)==daywind(2);
nsum1 = sum(L1);
i1=find(L1);
ibelow=i1-nwind+1;
if ibelow(1)<0;  % do not have full search window for first "year"
   i1(1)=[];
   ibelow(1)=[];
   nyr=nsum1-1;
else;
   nyr=nsum1;
end;
ihi=i1;
ilow=ibelow;
clear L1 nsum1 i1 ibelow;


% Check that start and end days for each year's search window are same
% month, dayofmonth, and julian day
i1 = X(ilow(1),[2 3 4]);
i2 = X(ihi(1),[2 3 4]);
I1=repmat(i1,nyr,1);
I2=repmat(i2,nyr,1);
J1=X(ilow,[2 3 4]);
J2=X(ihi,[2 3 4]);
if ~all(all(I1==J1)) | ~all(all(I2==J2));
   error('All start or end days of search windows are not same');
end;

% Store strings on search window and length of averaging period
str1 = sprintf('Search Window = %s %2.0f to %s %2.0f',...
   nmsmon{i1(1)}, i1(2), nmsmon{i2(1)}, i2(2)); 
str2 = sprintf('Averaging Period = %3.0f days',n);



% Loop over years
yr = repmat(NaN,nyr,1);
q  = repmat(NaN,nyr,1);
T  = repmat(NaN,nyr,4);

for m = 1:nyr;
   yr(m)= X(ihi(m),1); % year, corresponding to calendar year of end of search window
   iget = (ilow(m):ihi(m))'; % rows of X to pul 
   Y = X(iget,:); % Pull block of the tsm with the window for a year
   y = Y(:,5); %  daily data values
   [w,iw] = ndayma1(y,n); % returns n-day averages  w ending in sequential years iw
   %   indexed relative to rows of y, Y and X;
   % Compute low
   [z,iz]=nanmin(w); % the lowest n-day mean , and the index of the ending day of the period
   if isempty(iz);
       T(m,:)=[yr(m) NaN NaN NaN];;
   else;
       ix= iw(iz);  % row of Y holding ending day of n-day period with low flow
       q(m) = z;
       T(m,:)=Y(ix,[1 2 3 4]);
   end;
   
end;



