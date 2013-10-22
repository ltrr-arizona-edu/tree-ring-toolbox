function [A,B]=xdsub1(xend,yend,nwind,nlap,nmax)

% Form subscript matrices to be used in xdate.m.  
% Say you have two cores, a master of length 500 yr and another
% of length 250 yr.  Say you guess that the series are properly
% aligned when ring 500 of the master matches ring 241 of the
% other.  Calling function would cull out the two subseries,
% 1-500 of the master and 1-241 of the other and correlate 
% various overlapping segments at different lags.  This function
% makes arrays of subscripts to pull out correct segments of each
% series.

% D. Meko 7-14-93

% Use:  called by functions that attempt to aid in crossdating

% Other User-written functions needed:
%  None


%**********  INPUT ARGS
%
% xend -- value for this year (ubscript) of master series will be lined
%      up with value specified by yend for the undated sample series.
% yend -- corresp. subscript for y
% nwind -- width of window (number of years) for correlation analysis
% nlap -- overlap of windows (e.g., 10 years)
% nmax -- want no more than this number of windows returned.  The last
%     nmax windows will be given in A,B


%***************  OUTPUT ARGS
%
% A [nwind x mA] subscrpts of x giving years for correlation 
% B [nwind x mB] same for y



%***********  NOTES
%
% mA = mB
% 
% nlap can be as low as 1 to give maximum resolution in period of 
%  interest.
% xend and yend specify the alignment of the master and undated series.
%  The method used here is to specify the elements of x and y for the
%  alignment at the right-hand (most recent) side.



%*********************8   CODE


x=(1:xend)';   y=(1:yend)';   % cvs of sequential numbers

Lx=x-nwind >=0 ;  % only these elements could be valid ending seq nos
%		of nwind-year sub periods in x
ix=x(Lx);  % cull out these sequence numbers
ix= max(ix):-nlap:min(ix) ; % pull out a subset of sequence numbers of 
%	ending years, depending on the overlap specified
nx=length(ix);
ix = ix(nx:-1:1);  % make ix go from low to high

Ly=y-nwind >=0 ;  % only these elements could be valid ending seq nos
%		of nwind-year sub periods in y
iy=y(Ly);  % cull out these sequence numbers
iy= max(iy):-nlap:min(iy) ; % pull out a subset of sequence numbers of 
%	ending years, depending on the overlap specified
ny=length(iy);
iy = iy(ny:-1:1);  % make iy go from low to high



LL2= min([length(ix) length(iy) nmax]);  %  how many periods to get

ixx=ix(length(ix)-LL2+1:length(ix));  % ending subscript for each period
iyy=iy(length(iy)-LL2+1:length(iy));


% Dupe rows of rvs ixx, iyy.  Make a cv of length nwind.  Dupe cols of the cv.
% Then subtract arrays to get desired subscript arrays.
ixx=ixx(ones(nwind,1),:);
iyy=iyy(ones(nwind,1),:);
cv1=(nwind-1:-1:0)';
cv1=cv1(:,ones(LL2,1));
A=ixx-cv1;
B=iyy-cv1; 
