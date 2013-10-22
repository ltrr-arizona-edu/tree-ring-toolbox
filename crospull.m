function IX=crospull(nyears,nlag,plag)
% crosspull: make logical pointer to estimation and crossvalidation subsets for reconstruction model
% IX=crospull(nyears,nlag,plag);
% Last revised 7-13-99
%
% For a data matrix for climatic reconstruction with a model that possibly includes 
% lags on predictors, makes a logical pointer matrix to subsets of rows (years) to 
% be used in estimation and crossvalidation.  Uses "leave-n-out" strategy, where n 
% is two times the maximum of number of negative and positive lags for the reconstruction
% model.  
%
%*** INPUT ARGUMENTS
%
%  nyears (1x1)i  number of years in the calibration/verification period
%  nlag (1x1)i    number of negative lags on the predictors
%  plag (1x1)i    number of positive lags on the predictors
%
%
%*** OUTPUT ARGUMENTS
%
%  IX (m1 x n1)L   logical array, elements of a col refer to rows of data
%		matrix in the calling program. First column of IX has 1's for rows to 
%     be used 
%
%*** REFERENCES
%
% Crossvalidation is discussed by Michaelsen, J., 1987, Cross-validation in statistical 
% climate forecast models, J. of Climate and Applied Meterorology 26, 1589-1600.
%
% The "leave-n-out" strategy is used in:
% Meko, D.M., 1997, Dendroclimatic reconstruction with time varying subsets of tree 
% indices: Journal of Climate, v. 10, p. 687-696.
%
%*** UW FUNCTIONS CALLED -- none
%*** TOOLBOXES NEEDE -- none
%
%*** NOTES
%
% Makes an intermediate array of indexes to be used in cross
% validation in regression models.  Assume you are doing a regression
% analysis with crossvalidation by a "leave n out" strategy.  There
% are n years of data, and you will repetitively estimate a regression
% equation, leaving out enough years so that the predictand for the
% calibration period is independent of tree-ring data that might
% be related to the key year of the left-out sequence.
%
% Given the number of years of data you have, and the number of positive
% and negative lags on the predictor(s) in the model, this function builds
% a 1/0 logical pointer to your data telling which years to use in the 
% current model, and a subscript array telling which year is to be 
% estimated using the current model.
%
% Calling program would use column one of IX for first model:  rows with 1 would
% be used to fit the model, rows would zero would not be used for the fit, and
% the crossvalidation estimate is for year 1 of the nyears.  Next the calling program  
% would use column 2 of IX:  rows with 1 would be used for the fit, rows 
% with zero would not, and the crossvalidation estimate is for year 2. And so on 
% through all nyears of the data.
%
% It is necessary to insure that the data 
% for the central year of the "left-out" segment is independent of the 
% data used to estimate the parameters.  crospull.m handles this problem
% conservatively by assuming that 2 times the larger of plag and nlag  
% should be the number of years buffer on either side of the central year.
% Thus a model with 1 positive lag and 2 negative lags would look like a 
% model with 2 positive lags and two negative lags.
%
%
%  m1 equals nyears
%  n1 is the number of independent predictions that can be done in cross
%	validation.  n1=m1


maxlag=max(nlag,plag);
nzero = 1 + 4*(maxlag);  %  This number of years will be left out
%	of each regression (except a few) in the crossvalidation procedure.
%  The central
%	year of the left-out set will be the year predicted and checked
%	in crossvalidation.  The exceptions are the first and last 
%   2*maxlag data points, for which fewer than nzero need to be zeroed
%   out.

nmods=nyears;  % Number of regression models to estimate
n1= nmods;

IX = ones(nyears,nmods);  % Initialize logical array as all ones


%*****  Fill an intermediate subscript array

L1=(1:nzero)';  % cv like [1 2 3 4 5]'
L2=L1(:,ones(nmods-4*maxlag,1));  % dupe that column multiple times
K1= 0:(nmods-4*maxlag-1);  %  rv like [0 1 2 3 ...  45]
K2=  K1(ones(nzero,1),:);  % dupe that row multiple times
LK=L2+K2;  % For each regression model, a col of LK tells which rows are
%	to be left out.
%*******  Loop for "inside" years 

first=2*maxlag+1;
last=nyears-(2*maxlag);
for j=1:(nyears-4*maxlag);
	ithis=LK(:,j);
	IX(ithis,j+2*maxlag)=(zeros(nzero,1));
end



%***********   HANDLE "OUTER" YEARS

ncols=2*maxlag;  % this many columns are "outer" on each end
if ncols~=0;  % no outer-year problems if no lags
   xl=ones(nyears,ncols);
   for i=1:ncols;
      for j=1:(2*maxlag+i);
         xl(j,i)=0;
      end
   end
   xr=fliplr(xl);
   xr=flipud(xr);
   IX(:,1:(2*maxlag))=xl;
   IX(:,(nyears-2*maxlag+1):nyears)=xr;
end;

IX=logical(IX);
