function [est,r,N,xkey,rat2]=medrat(Y,X,overlap,yr,m)
% medrat: median ration estimation of missing monthly precipitation data
% CALL: [est,r,N,xkey,rat2]=medrat(Y,X,overlap,yr,m);
%
% Meko 1992, 1998
%
%****   INPUT ARGUMENTS  ******
%
% Y   (m1 x 13) monthly ppt array for station needing estimation
%		First col is year.
% X   (m2 x 13) monthly ppt array for station to be used as estimator
% overlap (1 x 2)  first and last year of overlapping period to
%		be used for computing median ratio.
% yr	year of value to be estimated
% m   month of value to be estimated
% xkey (1 x 1)r  predictor series value for the missing obs
% rat2 (1 x 1)r  ratio
%
%
%*******  OUTPUT ARGUMENT *******
%
% est  (1 x 1) the estimated value for year yr and month mo
% r (1 x 1)r  correlation coeff between the series for the 2 stns, computed
%   on the subset of years used to form the median ratio
% N (1 x 1)i  sample size (num of yr) for the median ration computation
%
%************ NOTES
%
% Source: Bradley (1976, p. 28)

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
%*****  missing data in Y or X;  and have zero values in some X.
%*****  Such years will be deleted before computing median ratio.
% 5-98: revised so that missing values must be NaN


LY1=isnan(Y1);
LX1=isnan(X1)  | X1==0;
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

rat1= Y1 ./ X1;
rat2 = median(rat1);
est = rat2 * xkey;

