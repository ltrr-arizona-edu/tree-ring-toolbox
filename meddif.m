function [est,r,N,xkey,dif2]=meddif(Y,X,overlap,yr,m)
% medrat: median difference estimation of missing monthly temperature data
% CALL: [est,r,N,xkey,dif2]=meddif(Y,X,overlap,yr,m);
%
% Meko  1998
%
%****   INPUT ARGUMENTS  ******
%
% Y   (m1 x 13) monthly tmp array for station needing estimation
%		First col is year.
% X   (m2 x 13) monthly tmp array for station to be used as estimator
% overlap (1 x 2)  first and last year of overlapping period to
%		be used for computing median difference
% yr	year of value to be estimated
% m   month of value to be estimated
%
%*******  OUTPUT ARGUMENT *******
%
% est  (1 x 1) the estimated value for year yr and month mo
% r (1 x 1)r  correlation coeff between the series for the 2 stns, computed
%   on the subset of years used to form the median difference
% N (1 x 1)i  sample size (num of yr) for the median difference computation
% xkey (1 x 1)r  tmp data for predictor station 
% dif1 (1 x 1)r  median difference, predictand station minus predictor
%
%************ NOTES
%
% Paradigm is median ratio method (Bradley 1976, p. 28,  and medrat.m)
%
% Idea is that the difference in temperature at two stations is conservative.
% So if temperature at station B is "typically" 5 degrees warmer than at station A, 
% a reasonable estimate of a missing monthly value at B is 5 degrees plus the value
% at A.  Median is used to identify "typical" difference.


I=find(X(:,1)==yr);  % index to year for estimation


% Form 0-1 vectors pointing to overlap period in Y,X 
LY = Y(:,1) >= overlap(1)  & Y(:,1) <= overlap(2);
LX = X(:,1) >= overlap(1)  & X(:,1) <= overlap(2);

% Pull out subset of overlap years for the specific month of est
% Compute ratio for each year; Compute median ratio.
% Mult X-value times ratio to get estimated value for Y

Y1 = Y(LY,m+1);
X1 = X(LX,m+1);

%****  New code (2-22-92) to allow period for long-term means to have 
%*****  missing data in Y or X; These years will not deleted from 
% subset of data used to compute median difference

LY1=isnan(Y1);
LX1=isnan(X1);
LL=LY1 | LX1;
Y1(LL)=[];
X1(LL)=[];

% correlation coef
r = corrcoef(X1,Y1);
r = r(1,2);

% Sample size
N = length(X1);

%**** end  new code
xkey = X(I,m+1);
dif1= Y1 - X1;
dif2= median(dif1);
est = dif2 + xkey;

