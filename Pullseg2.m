function [yr,Y]= pullseg2(x,yrx,len,offset)
% pullseg1: sequential len-year batches of time series x, with allowance for missing values
% [yr,Y]= pullseg1(x,yrx,len,offset);
% Last revised 1-1-01
%
% Low level function used prior to hinge.m by sprdloc.m 
%
%*** INPUT ARGS
%
% x (mx x 1) time series
% yrx (mx x 1) corresponding years for the time series x
% len (1 x 1) desired length of the segments (# of years)
% offset (1 x 1) each segment offset this many years from the previous
%
%
%*** OUTPUT ARGS *********************************
%
% yr (1 x ?) ending year of each len-year segment
% Y (mY x nY) each col is a time series segment, nY total segments
%     mY and nY computed in function
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS NEEDED --- none
%
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES *************************************
%
% Template was pullseg1.m, which does not handle missing values

% Remove any NaNs
L=isnan(x);
if any(L);
    x(L)=[];
    yrx(L)=[];
end;


% Preallocate and size.
mx = length(x);  % number of observations in input time series, after  removing NaNs

% Compute number of possible segments of the specified length and overlap
nseg1 = mx - len +1;
rs = mx:-offset:len;  % row subscript of last years of segments within NaN-cleaned x
yr = (yrx(rs))'; % rv of last year of each segment; segments not necessarily offset by same number of years
nseg2 = length(rs);  % number of possible segments 

g = [len-1:-1:0]';  % cv to be duped
G = repmat(g,1,nseg2);  % dupe cv g to matrix
%G = g(:,ones(nseg2,1));  % dupe cv g
RS = repmat(rs,len,1); % dupe rv rs to matrix
%RS = rs(ones(len,1),:);  % dupe rv rs

I = RS - G;  % index to rows of x

%yr = rs + yrx(1) - 1;  % rv of ending years of segments
X = x(:,ones(nseg2,1));  % dupe time series cv into matrix
Y = X(I); % pull out desired rows into  cols of Y
yr=fliplr(yr);  % make leftmost ending year the earliest segment
Y=fliplr(Y);  % likewise for the cols of Y