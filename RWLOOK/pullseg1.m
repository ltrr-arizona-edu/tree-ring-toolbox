function [yr,Y]= pullseg1(x,yrx,len,offset)

% len-year segments of time series x, offset len years
%
% D. Meko 12-13-93
%
%
%****************************  INPUT ARGS
%
% x (mx x 1) time series
% yrx (mx x 1) corresponding years for the time series x
% len (1 x 1) desired length of the segments (# of years)
% offset (1 x 1) each segment offset this many years from the previous
%
%
%********** OUTPUT ARGS *********************************
%
% yr (1 x ?) ending year of each len-year segment
% Y (mY x nY) each col is a time series segment, nY total segments
%     mY and nY computed in function
%
%
%*******  USER-WRITTEN FUNCTIONS NEEDED --- none
%
%
%******* GLOBALS   -- none
%
%
%******* NOTES *************************************
%
% A low level function to be used with rwlook and other functions. First
% written for rwlook to set up data for computation of Kendall's tau on
% different segments of time series of transformed ring width.



% Preallocate and size.
mx = length(x);  % number of years in input time series


% Compute number of possible segments of the specified length and overlap
nseg1 = mx - len +1;
rs = mx:-offset:len;  % row subscript of last years of segments
nseg2 = length(rs);  % number of possible segments 


g = [len-1:-1:0]';  % cv to be duped
G = g(:,ones(nseg2,1));  % dupe cv g
RS = rs(ones(len,1),:);  % dupe rv rs

I = RS - G;  % index to rows of x

yr = rs + yrx(1) - 1;  % rv of ending years of segments
X = x(:,ones(nseg2,1));  % dupe time series cv into matrix
Y = X(I); % pull out desired rows into  cols of Y
yr=fliplr(yr);  % make leftmost ending year the earliest segment
Y=fliplr(Y);  % likewise for the cols of Y