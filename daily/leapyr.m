function L=leapyr(yr)
% leapyr:  find out whether a year is a leap year (has a Feb 29)
% L=leapyr(yr);
% Last revised 2-14-01
%
% Operates on scalar or column vector of years.  Returns whether each year is a leap year.
%
%*** INPUT
%
% yr (? x 1)i  year or years to be tested; e.g., [1800; 1900 ; 2000];
%
%
%*** OUTPUT
%
% L (? x 1)L  1 if year in yr is a leap year, 0 otherwise
%
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% Feb of a leap year has 29 days.  Feb of a non-leap year has 28 days.  Leap year is every 
% year evenly divisible by 4, with special condition on century years (e.g., 1900, 2000).
% If century year not evenly divisible by 400, NOT a leap year.

L=rem(yr,4)==0  & ~(rem(yr,100)==0 &  rem(yr,400)~=0);