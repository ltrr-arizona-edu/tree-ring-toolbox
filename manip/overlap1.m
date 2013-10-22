function [Y,yrY,yrnan]=overlap1(x1,yrx1,x2,yrx2,kopt)
% overlap1:  overlapping data between two time series
% Last revised 5-17-00
% [Y,yrY]=overlap1(x1,yrx1,x2,yrx2,kopt)
%
% Overlapping data between two time series. Given two time series and their year 
% vectors, returns the 2-column matrix of data for years in common, and the 
% year vector for the matrix.  Optionally can include only non-NaN data in the matrix.
%
%*** INPUT ARGUMENTS
%
% x1 (mx1 x 1)r   time series 1
% yrx1 (mx1 x 1)i   year vector for x1
% x2 (mx2 x 1)r time series 2
% yrx2(mx2 x 1)i year vector for x2
% kopt (1 x 1)i  options
%   kopt(1) ....   control over NaNs
%		==1 no special treatment of NaNs
%		==2 observations with NaN in either series omitted from returned data
%
%*** OUTPUT ARGUMENTS
% 
% Y (mY x 2)r   overlap data for x1 (col 1) and x2 (col 2)
% yrY (mY x 1)  year vector for Y (see notes)
% yrnan(? x 1)  years not in Y and yrY because data NaN at either site
%		[] if kopt(1)==1  or if kopt(1)==2 but no NaNs in the overlap
%
%*** REFERENCES -- NONE
%*** TOOLBOXES NEEDED -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%
%*** NOTES
%
% yrY.   yrY  and Y contain all overlapping years, except that observations with NaN
% in either x1 or x2 are omitted under control of kopt(1);
%
% Minimum overlap.  Error if fewer than 2 yr.


% Check inputs
[mx1,nx1]=size(x1);
[mx2,nx2]=size(x2);
if nx1 ~= 1 | nx1~=1;
   error('x1 and x2 should be cv');
end;
if mx1 ~= length(yrx1) | mx2 ~= length(yrx2);
   error('yrx1 and yrx2 should be same row size as x1, x2');
end;

yrgo=max([min(yrx1)   min(yrx2)]);
yrsp=min([max(yrx1)   max(yrx2)]);
if yrsp-yrgo<1;
   error('Less than 2 years of overlap');
end;

% Pull common data
L1=yrx1>=yrgo & yrx1<=yrsp;
L2=yrx2>=yrgo & yrx2<=yrsp;
Y=[x1(L1)  x2(L2)];

% Set common year vector
yrY=(yrgo:yrsp)';
yrnan =[]; % initialize the NaN-year vector

% Control for NaN
L3=isnan(Y);
L4=(any(L3'))'; % cv marking year with NaN in either
if any(L4);
   yrnan=yrY(L4);
end;


% Removal of NaN years
if kopt(1)==2 & ~isempty(yrnan);
   yrY(L4)=[];
   Y(L4,:)=[];
   if size(Y,1)<2;
      error('Y for common period after deleting NaN years has fewer than 2 rows');
   end;
end;
      
   





