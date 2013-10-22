function [cant,L1]=daymiss(X,S,xmiss)

% Is data coverage by predictor stations sufficient to fill in the missing
% daily precip values?
%
% D Meko 2-10-94
%
%************  INPUT ARGS
%
% X (mX x nX) daily precip matrix.  Cols 1,2 give year, day, 
%      Number of stations is nX-2.  mX equals 366 times the number of
%      years of data in X.
%
% S (mS x nS) index to stations in X to be used as predictors for
%		each station in X.  Each row for a key station.  Each col
%		for one of nS predictor stations. For example,
%		S(1,:)=[2 3 4] means stations 2,3, and 4 (cols 4,5,6) of
%		X are the predictor series
% xmiss (1 x 1) missing value code -- like -1.111
%
%
%********  OUTPUT ARGS
%
% cant(nS x 1) how many values can't be estimated for each station.
%		The ideal result is that cant is a column of zeros.
% L1(mX x nX) which rows of X hold the values that cannot be extimated
%		Col 1 is for first station, col 2 for second station, etc
%
%
%*******************  USE
%
%	Before dayfill.m to check that the predictor series you have in mind
%  have adequate data coverage to fill in all missing values in each
%  key series.  Use as follows:
%
%  Make sure X, xmiss, S are in workspace
%  Call daymiss.m
%  Check col vector "cant" for any nonzero values.  A 5 in 
%      row 3, say, means that series 3 has five missing values that
%		cannot be handled
%	Display a matrix of years, days with the problem values in series k  by
%		typing 
%
%		[X(L1(:,k),1)  X(L1(:,k),2)]


yr=X(:,1);
day=X(:,2);

X(:,1:2)=[];  % cut off year and day columns
[mX,nX]=size(X);
cant=zeros(nX,1);
L1=zeros(mX,nX);  % will hold 0/1 to rows unable to estimate
[mS,nS]=size(S);

L=X==xmiss;  % logical -- 1 if  a missing value

for i=1:nX;  % loop over all stations
	s=S(i,:);   % pointer to cols of predictor stations in X
	Lk=L(:,i);  % value missing at key
	Lp=L(:,s);  % value missing at predictors

	Lpall=  (all(Lp'))';  % cv, 1 if missing at all predictors stations
		% 0 otherwise

% Missing at key and all predictor stations, and not Feb 29 of a leap year
	L2=Lk & Lpall & day~=60 & rem(yr,4)~=0;  
	L1(:,i)=L2;
	cant(i)=sum(L2);  % sum of values missing at key and at all predictor stations
end



