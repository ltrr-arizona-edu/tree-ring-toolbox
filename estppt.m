function [Z,U,V,r,rsize]=estppt(Y,X,medpd,W,crit)

% Estimation of multiple months of missing monthly ppt value by median-
% ratio method (Bradley 1976, p. 28)

%****   INPUT ARGUMENTS  ******
%
% Y   (m1 x 13) monthly ppt array for station needing estimation
%		First col is year.
% X   (m2 x 13) monthly ppt array for station to be used as estimator
% medpd (1 x 2)  first and last year of period to
%		be used for computing median ratio.
% W (? x 2) month, year of data needing estimation .  First 2 rows
% 	have special information (see findmiss.m)
% crit  (1 x 1) missing value code numeric

%*******  OUTPUT ARGUMENT *******
%
% Z (m1 x 13)  Y, except with estimates replacing missing value codes
% U (? x 5) Record-keeping array telling which predictor series, 
% 	median-ratio period, etc were used.  If estimate other missing
%	values later from another station or using another medians 
%	period, can append to U.  Each row of U corresponds to a value
% 	of Y successfully estimated from X.  The cols are as follows:
%
% V (? x 2) information array like input array W, except including only
%   those months, years that could not be predicted from the current run.
%   For example, the predictor may also have missing values in the months
%   to be estimated.  The array V can then serve as the starting point
%   by acting as W in another call to estppt using a different predictor
%   series. 
%
%	1,2	month, year for an estimate
%	3	numeric id for predictor station
%	4,5	period used for computing median ratios
%
% r (1 x 12)  median ratios for months 1 (jan) thru 12 (dec) for the 
%   	      period meanpd
% rsize (1 x12) number of years of ratios available for each month
%	for forming the median-ratios in r.  This number can
%	differ from the number of years specified by medpd because
%	(1) zero-ppt years in X are unusable -- would give infinite ratio
%	(2) some months needing estimates may also have missing data in X


%*********************   PREALLOCATE AND INITIALIZE

clear V;
U=zeros(302,5);  % Will hold listing of which months estimated, and
         	 % from which predictor series.  Row size of 302
		 % allows for max of 300 estimated values for a given
		 % station.



Z=zeros(Y); % initialize what will be the estimated version of Y

idpred=input('NUMERIC ID NUMBER FOR PREDICTOR SERIES:    ');
%	This id number will be put in col 3 of U as a record for future
%	reference on how values estimated.

[mx,nx]=size(X);  % size of predictor series


% Form 0-1 vectors pointing to period for medians computation in Y,X.

LY = Y(:,1) >= medpd(1)  & Y(:,1) <= medpd(2);
LX = X(:,1) >= medpd(1)  & X(:,1) <= medpd(2);

% Pull out subset of medpd years for the specific month of est
% Compute ratio for each year; Compute median ratio.

for  m = 1:12;  %  compute median ratios for each month

	Y1 = Y(LY,m+1);
	X1 = X(LX,m+1);

	LY1 = Y1 >= crit ;
	LX1 = X1 >= crit | X1 <=0;  
	LL = (LY1 | LX1);  %  ones if either X or Y missing data,
				% or X has zero value
	
	Y1(LL)=[];
	X1(LL)=[];

   rat = Y1 ./ X1;
   r(m) = median(rat);  % median ratio, month m
   rsize(m) = length(Y1);
end


[m2,n2]=size(W);  % input mtx of mos, yrs mising
WW=W(3:m2,:);   %Because first two rows have header info
m3=m2-2;  % number of values needing estimation

j=0;   % Initial index for U
V(1:2,1:2)=W(1:2,1:2);  % Array pointing to values that
	% cannot be estimated.

for i=1:m3;  %  Loop for each monthly value to be estimated
	month=WW(i,1); % Month of year for the estimate (1=jan)
	year = WW(i,2); % Year with value to be estimated
	rat2 = r(month);   % appropriate ratio for this month of year
	iiy = find(Y(:,1)==year);  % appropriate row subscript 
		%  of Y for this year
	iix= iiy + (Y(1,1) - X(1,1)); % corresp row subscript for Y
	yvalue=Y(iiy,month+1);  % Val of predictand -- should be missing


	if iix <1 | iix > mx;  % Predictor data period does not cover this year!
		est=yvalue;  % X data value cannot be estimated.
		V=[V ; month  year ] ;

	else;  % Predictor series covers this year
		xvalue=X(iix,month+1);  % Value predictor for the missing month

	   if xvalue >= crit;  %  If X has missing data in this month, Y
		% missing value cannot be replaced. Set Y value as
		% its old value (missing).

		V=[V ; month  year ] ;
		est=yvalue
	   else
   		est = rat2 * xvalue;
		j=j+1;
		U(j+2,1:2)=WW(i,1:2);
	   end
	end

	Z(iiy,month+1) = est;
	
end

% Make the array recording the estimation info for this run

U=U(1:j+2,:);  % cut off unused rows of initialized U
U(3:j+2,3:5)=[idpred(ones(j,1),:)   medpd(ones(j,1),:)];
U(1:2,1:2)=W(1:2,1:2);




 
