function [Y,N,yrY,R]=lfsync(X,I1,b,yrs,xmiss,k)

% Summarize synchrony in low-frequency time series variations
%
% D. Meko, 2-21-94
%
%
%**********  INPUT ***************************
%
% X (mX x nX) matrix of time series, end parts may have missing data
%    Col 1 holds year.  Col 2-nX holds time series for nX-1 variables.
%   Each row a year.  mX years might be larger than desired analysis period.
% I1 (mI1 x 1)    Index to cols of X specifying which cols to use in analysis.
%   Year column skipped in indexing:  e.g., I1(1)=4 means col 5 of X.
% b (1 x nb)  lowpass filter weights; nb must be odd.
% yrs (1 x 2)  period for analysis;  must be included within X.
% xmiss (1 x 1) missing-value code (e.g., 99)
%
%
%***********************  OUTPUT  ****************************
%
% Y (mY x 5) time series of prob points for each year of smoothed key period.
%   mY rows corresponding to years;  The mY years will be a subset of 
%   period "yrs", truncated because lose ends in smoothing.
%
%	  col 1 = minimum
%    col 2,3,4 = .25, .50, .75 quantiles
%    col 5 = maximum
%
% nY (mY x 1) sample size in each year of Y.  Sample size is number of
%   trees (or cores) in each year.
% yrY (mY x 1) year vector for Y
%
%
%********************** USER-WRITTEN FUNCTIONS NEEDED **********  NONE
%
%
%***************  USE  *******************************************
%
% Three tree-ring sites in northern Great Plains.  Culled out very long core or
% tree indices for subset of 15 series.  Observed in fitting growth trends that
% some very long wavelength variations.  Correlation of indices with climate
% in this century unlikely to show whether the low-freq variations real or not.
% So check for internal consistency using nonparametric (rank) methods.
%
%
%**************  STEPS  ***********************************
%
%  *	Cull out subset of X for analysis
%  *	Smooth the series, which will also lop off some data from ends
%  *	Loop for each series
%    *	Compute rank within nonmissing years in key period for each 
%		year of data
%  *	Loop for each year in key period
%    *  Get sample size, extremes, 25, .5, .75 quantiles for each year
%*******************************************************


%*****  SIZE AND PREALLOCATE

if nargin~=6, error('nargin not equal to 5'), end
[mX,nX]=size(X);
if yrs(1)<X(1,1), error('yrs(1) before first year in X'), end;
if yrs(2)>X(mX,1), error('yrs(2) after last year in X'), end;
[mI1,nI1]=size(I1);

[mb,nb]=size(b);
if mb~=1 & rem(nb,2)~=1, error('b not an odd-length row vector'), end

if nI1~=1, error('I1 must be col vector'), end;
if max(I1) > nX-1, error('I1 points to too high a column for X'), end;
if min(I1) < 1, error('I1 must point to col 2 or higher of X'), end;


%***** Cull out subset of X for analysis


L1=X(:,1)>=yrs(1) & X(:,1)<=yrs(2);
X=X(L1,[1 (I1+1)']);  % +1 skips year column 
[mmX,nnX]=size(X);


% allocate space for rank matrix R
if k==1
	R=xmiss(ones(mmX-nb+1,1),ones(mI1+1,1));   % initially fill with missing values
	yrY=((yrs(1)+fix(nb/2)):(yrs(2)-fix(nb/2)))';
elseif k==2
	R=xmiss(ones(mmX,1),ones(mI1+1,1));   % initially fill with missing values
	yrY=(yrs(1):yrs(2))';
end


R(:,1)=yrY;


% Smooth and form nonexceed prob time series of each variable


for i=1:mI1;  % loop for every series
	z=[];
	x1=X(:,i+1);
	L2=x1~=xmiss;  % logical pointer to nonmissing data rows
	x2=x1(L2);
%	[y2,yry2]=mafilt2(x2,X(L2,1));  % smooth the series
	
	[y2,yry2]=filter1(x2,X(L2,1),b,2);

	L3=yrY>=yry2(1) & yrY<=yry2(length(y2));  % for replacing in R
	[y2s,I3]=sort(y2);
	z(I3)=(1:length(y2))'/(length(y2)+1);


	R(L3,i+1)=z';
end


% Lop off beginning and ending years that have no quantiles
[mR,nR]=size(R);
RR=R(:,2:nR);
RR=RR';
RL=all(RR==xmiss);
RL=RL';
R(RL,:)=[];
[mR,nR]=size(R);

Y=xmiss(ones(mR,1),ones(6,1));
Y(:,1)=R(:,1);


%  Fill the quantiles matrix by looping through by year

[mY,nY]=size(Y);
N=zeros(mY,1);
for i=1:mY;
	disp(i);
	yy=R(i,2:nR);
	Ly=yy~=99;
	Lys=sum(Ly);
	N(i)=Lys;
	if Lys==0
		error('all 99 in this year')
	else
		v=yy(Ly);  % the sample
	end;
	if Lys==1;  % can only get median
		Y(i,4)=v;
	end
	if Lys==2 | Lys==3;  % can get median and extremes only
		Y(i,2)=min(v);
		Y(i,6)=max(v);
		Y(i,4)=median(v);
	end
	if Lys>=4;
		Y(i,2)=min(v);
		Y(i,3)=quantile(v',.25);
		Y(i,5)=quantile(v',.75);
		Y(i,6)=max(v);
		Y(i,4)=median(v);
	end

end
