function nd=yr2seq(day,mon,yr,yrbeg)
% Compute sequential number of day, assuming day 1 is Jan 1, yrbeg
%  where yrbeg is divisible by 4.


still not correct

%*********  INPUT ARGS
%
% day -- day of month
% mon -- month of year (1=jan)
% yr -- year (e.g., 1976)
% yrbeg -- arbitrary beginning (leap) year for sequence. 
%  Observation 1 is for January 1 of year yrbeg.


%*****************  OUTPUT ARGS
%
% nd -- the sequence number for the day


% Make sure beginning year is a leap year.
if(rem(yrbeg,4)~=0), stop('yrbeg is not divisible by 4'); end

% How many previous full years?
n=yr-yrbeg;

% How many days in those years?
nleap=fix(n/4)+1;  % number of leap years in those years
if n ~=0
	ndl=365*n+nleap;   % number of days in the lead-in period.
else
	ndl=0;  % If no lead in period
end


% What  is sequence number of day,mon,yr? 
% Start with cumulative sum of a non leap year.
cumday1=cumsum([0 31 28 31 30 31 30 31 31 30 31 30]);
cumday2=cumsum([0 31 29 31 30 31 30 31 31 30 31 30]);
if rem(yr,4) == 0;   % year is a leap year
	cumday=cumday2;
else
	cumday=cumday1;
end

nd=ndl+cumday(mon) + day;

