function jday=mdy2jdy(month,day)
% mdy2jdy: convert month, day, to julian day (sequential day in 366-day year)
% CALL: jday=mdy2jdy(month,day);
%
% Meko 5-12-98
%
%********************  IN 
%
% month (1 x 1)i month of year (1==jan)
% day (1 x 1)i  day of month 
%
%***************** OUT 
%
% jday (1 x 1)i  Julian day (Jan 1 ==1, Dec 31 == 366)
%
%
bef = [0 31 29 31 30 31 30 31 31 30 31 30];
bef = cumsum(bef); % increment of days 


jday = day + bef(month);
