function R = corr1(C,D,L1,endmo,k)

% Returns correlation coef matrix between seasonal sums or averages of 
% climate variables and a response variable -- like a tree-ring index

% Or returns the array of first order autocorrelation coeffs of the 
% precip or temp series

% Let's say the climate variable is monthly precip and the response 
% variable is the tree-ring index.

%****************  ARGUMENTS ********************************

% D - monthly ppt array. Each row a year t.  Each col a month.
% 	36 cols represent three years.   Col 1 is the year t.  D can cover
%  a longer period (number of rows) than the correlation analysis will
%  cover.
% 		col 2   jan t-2
% 		col 3	feb t-2
%		col 13	dec t-2
%		col 14  jan t-1
%		col 26  jan t
%      col 37  dec t

% C - the tree-ring series for analysis in this call (cv).
%	 either as-is or prewhitened,
% 	 and includes only the years specified in L2
% L1 - logical array specifying rows (years) of D for correlation analysis

% endmo - ending month of tree-ring growth year.  1=jan, 12=dec.
%	So if want to gear analysis to growth year ending in august of
%	year t, set endmo=8.

% k - type of analysis 
% 	1 ppt vs tree
%	2 temp vs tree
%	4 ar(1) coef ppt
%	5 ar(1) coef temp

%********************  CAUTIONS ******************************

% D should have been built so that cols 26-37 correspond to jan-dec of
%   year t.  
% D should not have missing values over the year range specified by
%   L1
% 
%********************  OUTPUT   *********************************

% R - the correlation matrix, 12 x  12, between tree rings and seasonal
%	ppt

% The ending month of the monthly grouping varies down the rows.
% The number of months in the grouping varies across the cols.

% Example:  assume endmo=8.

% Col 1 contains r between tree-rings and one-month total ppt
%	col 1, row 1 is for August of year t
%	col 1, row 2 is for July   of year t
%     col 1,row 12 is for Sept   of year t-1
% Col 2 contains r between tr and two-month total ppt
%	col 2, row 1 is for July-Aug, year t
%	col 2, row 2 is for June-July, year t
% 	col 2,row 12 is for Aug t-1 thu Sept t-1
% Col 12 contains r between tr and 12-month total ppt
%	col 12, row  1 is for Sept t-1 thru Aug t
%	col 12, row  2 is for Aug  t1  thru July t
%     col 12, row 12 is for Oct t-2 thru Sept t-1

%*********  PREALLOCATE *********************************************

R = zeros (12);  T=zeros(13,13);

J=[1:13:144]';

S5=zeros(sum(L1),12);

lastmo=37- (12-endmo);  % tie ending month to column of D

for k1=1:12;  % loop for each number of months in seasonal grouping
	top=lastmo-(k1-1):lastmo;  % rv.  Say [32 33] for k1 = 2,endmo=8
	S1 = top(ones(12,1),:); % dupe rows of top
	S2 = (0:11)';
	S3 = S2 (:,ones(k1,1)); % dupe cols of S2
	S4 = S1-S3;

	for k2 = 1:12;  % sum or ave the clim variable over months
		if k1==1; % special case for vector
			S5(:,k2) = D(L1,S4(k2,:));
		else 
			S5(:,k2) = (sum((D(L1,S4(k2,:)))'))';
		end ;  % of if loop
	end;  % of for loop k2
	S5=dtrend(S5);

	if k==1 | k==2; % trees vs ppt or trees vs temp
		T = corrcoef ([C  S5]);
		R (:,k1) = (T(1,2:13))';
	elseif k==4 | k==5; % ar(1) coefs of PPt or temp
		S6=covf(S5,2);
		S6=S6(J,:);
		R(:,k1)= S6(:,2) ./ S6(:,1);
	else
	end

disp(k1)

end; % of loop for k1
 
