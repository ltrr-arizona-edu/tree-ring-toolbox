function IX=crospul1(nyears,nlag,plag)

% D. Meko 8-20-95; revised form of crospull.m .  Revised to be
% less stringent in "leaving out" years around the validation
% target

% Makes an intermediate array of indexes to be used in cross
% validation in regression models.  Assume you are doing a regression
% analysis with crossvalidation by a "leave n out" strategy.  There
% are n years of data, and you will repetitively estimate a regression
% equation, leaving out enough years so that the predictand for the
% calibration period is independent of tree-ring data that might
% be related to the key year of the left-out sequence.

% Given the number of years of data you have, and the number of positive
% and negative lags on the predictor(s) in the model, this function builds
% a 1/0 logical pointer to your data telling which years to use in the 
% current model, and a subscript array telling which year is to be 
% estimated using the current model.

% It is necessary to insure that the data 
% for the central year of the "left-out" segment is independent of the 
% data used to estimate the parameters.  crospul2.m handles this problem
% by assuming that 1 times the larger of plag and nlag  
% should be the number of years buffer on either side of the central year.
% Thus a model with 1 positive lag and 2 negative lags would look like a 
% model with 2 positive lags and two negative lags.  A model with 1 
% pos and 1 neg lag would be a "leave-3-out" model

%*******   INPUT ARGUMENTS
%
%  nyears (1x1)  number of years in the calibration/verification period
%  nlag (1x1)    number of negative lags on the predictors
%  plag (1x1)    number of positive lags on the predictors


%*********  OUTPUT ARGUMENTS
%
%  IX (m1 x n1)   logical array, elements of a col refer to rows of some
%		matrix in the calling program


%********  NOTES
%
%  m1 equals nyears
%  n1 is the number of independent predictions that can be done in cross
%	validation.  n1=m1


maxlag=max(nlag,plag);
nzero = 1 + 2*(maxlag);  %  This number of years will be left out
%	of each regression (except a few) in the crossvalidation procedure.
%  The central year of the left-out set will be the year predicted and 
%  checked in crossvalidation.  The exceptions are the first and last 
%  maxlag data points, for which fewer than nzero need to be left
%  out.

nmods=nyears;  % Number of regression models to estimate
n1= nmods;

IX = ones(nyears,nmods);  % Initialize logical array as all ones


%*****  Fill an intermediate subscript array

L1=(1:nzero)';  % cv like [1 2 3 4 5]'
L2=L1(:,ones(nmods-2*maxlag,1));  % dupe that column multiple times
K1= 0:(nmods-2*maxlag-1);  %  rv like [0 1 2 3 ...  45]
K2=  K1(ones(nzero,1),:);  % dupe that row multiple times
LK=L2+K2;  % For each regression model, a col of LK tells which rows are
%	to be left out.
%*******  Loop for "inside" years 

first=maxlag+1;
last=nyears-(maxlag);
for j=1:(nyears-2*maxlag);
	ithis=LK(:,j);
	IX(ithis,j+1*maxlag)=(zeros(nzero,1));
end



%***********   HANDLE "OUTER" YEARS

ncols=maxlag;  % this many columns are "outer" on each end
if ncols~=0;  % no outer-year problems if no lags
xl=ones(nyears,ncols);
for i=1:ncols;
	for j=1:(maxlag+i);
		xl(j,i)=0;
	end
end
xr=fliplr(xl);
xr=flipud(xr);
IX(:,1:(maxlag))=xl;
IX(:,(nyears-maxlag+1):nyears)=xr;
end;
