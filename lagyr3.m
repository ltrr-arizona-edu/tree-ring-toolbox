function [Z,yrZ] = lagyr3(X,yrX,lg)
% lagyr3:  build a lagged predictor matrix for a time series matrix
% CALL: [Z,yrZ] = lagyr3(X,yrX,lg);
%
% D. Meko 11-18-92, revised 3-15-98
%
%%******   INPUT ARGS
%
% X (mX X nX) unlagged array of predictors, mX years, nX variables
% yrX (1 x 2)  beginning and ending years of X
% lg (3 x 1)  shifting and lagging parameters
%		1 number of years to shift X (the delay; always negative)
%		2 number of negative lags to include on X
%		3 number of positive lags to include on X
%
%*****   OUTPUT ARGS
%
% Z (mZ x nZ) lagged and possibly shifted predictor matrix
%		mZ=mX  same number of rows (years) as original matrix
%		nZ= (1+lg(2)+lg(3)) * nX    depends on lag parameters
% yrZ (3x2) 
%     row 1= "years" corresponding to rows of lagged array Z.  
%		Function designed such that yrZ would be years of reconstructed
%		variable were the matrix Z used in regression.  If no delay 
%		(i.e., if lg(1)=0),  yrZ is identically equal to yrX.  Depending
%		on values of lags, however, some leading or ending rows of Z
%		may contain NaNs.  
%     
%		row 2 = totally valid years - corresp to rows of Z with no NaNs
%			Calling pgm could use this information to automatically
%         Determine possible years for reconstruction if use full
%         number of specified lags and delays in regression model.
%
%		row 3 = row indices of Z corresponding to the years in row 2
%
%
%******************* NOTES **********************
%
% Assumes matrix of predictor variables, specified maximum number of
% positive and negative lags that might be needed on those variables, and
% optional delay in the response.
%
% Puts unlagged predictors in first cols of Z; lags -1, -2 ... in
% next cols; and lags +1, +2 .... in remaining cols.  
%
% Input X array covers specified years.  Output Z array is same row
% size, and can cover same period (if no delay) or the original period
% shifted back by some number of years equal to the specified delay.  
%
% For the most blind use of this function, pass the predictor array and 
% lag/shift parameters.  Get back the Z array of lagged/shifted
% predictors.  yrZ(2,:)  gives years you can reconstruct if using all
% lags and any specified delay.   yrZ(3,:) gives the corresponding rows
% of Z to pull off as the predictor data. 
%
% A less automatic use of the function is to not use all the possible
% lags (lg) in the regression model.  In this case, may be able to use more
% years than given by yrZ(:,2).  For example, maybe you already have
% called this function  with lg= [0 1 1], but later decide to build
% a regr model with only a negative-1 lag.  The calling program can figure
% out the appropriate usable rows of Z.  You can double check by making
% sure no NaNs in whatever subset array of Z you extract as predictor
% matrix.
%
% The array Z is intended for dual purpose:  as a mother array from which
% to pull subsets of rows and cols for the predictor matrix  in
% regression;  and as the mother array to pull the same subset of cols,
% but a bigger set of rows to substitute into the regression
% model to get the long-term reconstruction.


%**********   OTHER USER-WRITTEN FUNCTIONS NEEDED  -- NONE



%*********  CHECK NUMBER OF ARGUMENTS

L1=[nargin nargout];
if (nargin~=3 & nargout~=2), 
	error('REQUIRED:  NARGIN=3,  NARGOUT=2');
end

if(lg(1)>0 | lg(2)<0  | lg(3)<0),
	error('lg(1) must be non-positive, lg(2) and lg(3) non-negative')
end

%******   PREALLOCATE, SIZE, INITIALIZE

yrZ=zeros(3,2);

[mX,nX]=size(X);
if(mX~=(yrX(2)-yrX(1)+1)), 
	error('N of YEARS FOR INPUT ARRAY INCONSISTENT'), 
end

a=NaN;
mX1=mX+lg(2)+lg(3)-lg(1);  % number of rows in intermediate array that
%	has NaN rows possibly tacked onto beginning and ending rows of X
X1=a(ones(mX1,1),ones(nX,1));  % that intermediate array init to NaN

mZ=mX;  % output array Z same number of years (rows) as X
nZ=(1+lg(2)+lg(3))*nX;  % numb of cols in output array Z
Z=a(ones(mZ,1),ones(nZ,1));  % an output argument


%******  BUILD X1

I1=((lg(2)+1):(mX+lg(2)))';  % index to rows of X1 to be replaced with X
X1(I1,:)=X;


%*********  FILL Z WITH CORRECT SUBSETS OF X1

IZgo=lg(2)+1 ; %-lg(1);  % index for starting row of zero lag data
IZ=(IZgo:IZgo+mX-1);  % row indices to pull for zero-lag part
Z(1:mX,1:nX)=X1(IZ,:);

%  loop for negative lags
for i=1:lg(2)
	IN=IZ-i;
	Jgo=1+(nX*i);  % starting col of Z that this block goes into
	Jstop=Jgo+nX-1;
	Z(:,Jgo:Jstop)=X1(IN,:);
end

%  loop for positive lags
for i=1:lg(3)
	IP=IZ+i;
	Jgo=1+(nX*lg(2)) + (nX*i);
	Jstop=Jgo+nX-1;
	Z(:,Jgo:Jstop)=X1(IP,:);
end


%****************   WHAT YEARS DO ROWS OF Z MATCH UP TO?
%
% Normally, the delay (lg(1)) is zero, and the years of Z are
% considered to match the years of X.  But if lg(1) is non-zero,
% you want year to shift the year coverage of Z relative to some
% predictand series in the calling program.

yrZ(1,:)=yrX+lg(1);  % note that no effect if lg(1)=0

begvalid= yrX(1) + lg(1) + lg(2);
endvalid= yrX(2) + lg(1) - lg(3);

beg1 = 1 + (begvalid-yrZ(1));
end1 = beg1 + (endvalid-begvalid);

yrZ(2,:)=[begvalid endvalid];
yrZ(3,:)=[beg1 end1];
