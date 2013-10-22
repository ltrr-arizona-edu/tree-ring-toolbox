function [y,irow]=ndayma1(x,nday);
% ndayma1:  n-day moving average of a time series
% [y,irow]=ndayma1(x,nday);
% Last Revised 9-14-00
%
% For a time series x, computes the n-element (e.g., n-day) moving average and
% returns the result in y, along with a row index irow that is the ending 
% row indexed to x of the n-element average.
%
%*** INPUT
% 
% x (mx x 1)r   time series (e.g., daily precip);
% nday (1 x 1)i  number of days in averaging period; if nday==mx, a single value returned
%
%
%*** OUTPUT
%
% y (my x 1)r   time serie so n-day moving average of x, indexed to ending day
% irow (my x 1)r row (in x) corresponding to ending day of n-day average stored in y
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% NaNs are ignored in averages (function uses nanmean)
% nday:  If nday =mx, computes the simple average of the values in x and
%   returns this scalar as y along with scalar irow that is the last row-index of x


[mx,nx]=size(x);
if nx~=1;
   error('x must be cv ');
end;

% Special case: 1-day "moving average"
if nday==1;
   y=x;
   irow=(1:length(x))';
   return;
end;

% Special case: 366-day moving average -- this is just mean daily value for year
if nday==366;
   y = nanmean(x);
   irow = length(x);
   return;
end;

% Build row index to n-day periods
nper = mx-nday+1; % number of n-day periods 
j1 = (1:nday)'; % cv of row index into x for first n-day period
J1 = repmat(j1,1,nper);
jadd = [0:(nper-1)]; % rv increment
Jadd = repmat(jadd,nday,1); % dupe to matrix
J = J1+Jadd;  % matrix of row indices to x

% Pull matrix of n-day periods
X=x(J);

% Compute means
y = (nanmean(X))'; % n-day means, as cv
irow = (J(nday,:))'; % row index into x corresponding to ending day of n-day period



