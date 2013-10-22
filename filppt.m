%function [Z,r]=filppt(Y,X,meanpd,W,crit)

% Estimation of multiple months of missing monthly ppt value by median 
% ratio method
% Bradley (1976, p. 28)


%****   INPUT ARGUMENTS  ******
%
% Y   (m1 x 13) monthly ppt array for station needing estimation
%		First col is year.
% X   (m2 x 13) monthly ppt array for station to be used as estimator
% meanpd (1 x 2)  first and last year of period to
%		be used for computing median ratio.
% W (? x 2) month, year of data needing estimation .  First 2 rows
% 	have special information (see findmiss.m)
% crit  (1 x 1) missing value code numeric

%*******  OUTPUT ARGUMENT *******
%
% Z (m1 x 13)  Y, except with estimates replacing missing value codes
% r (12 x 1)  median ratios for months 1 (jan) thru 12 (dec) for the 
%   period meanpd
%

Z=Y; % initialize what will be the filled in version o f Y

% Form 0-1 vectors pointing to period for medians computation in Y,X 

LY = Y(:,1) >= meanpd(1)  & Y(:,1) <= meanpd(2);
LX = X(:,1) >= meanpd(1)  & X(:,1) <= meanpd(2);

% Pull out subset of meanpd years for the specific month of est
% Compute ratio for each year; Compute median ratio.
% Mult X-value times ratio to get estimated value for Y

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

end


 	[m2,n2]=size(W);  % input mtx of mos, yrs mising
 	WW=W(3:m2,:);   %Because first two rows have header info
	m3=m2-2;  % number of values needing estimation

for i=1:m3
	rat2 = r(WW(i,1));   % appropriate ratio for this month of year
	ii = find(X(:,1)==WW(i,2));  % appropriate row subscript for the year
   est = rat2 * X(ii,m+1);
	Z(ii,m+1) = est;
end


 
