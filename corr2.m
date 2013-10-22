function [xM,xS] = corr2(D,E,L1,endmo)

% Compute means and standard devs of various monthly groupings of 
% ppt or temp


%****************  ARGUMENTS ********************************

% D - monthly ppt array. Each row a year t.  Each col a month.
% 	36 cols represent three years.   Col 1 is the year t.
% 		col 2   jan t-2
% 		col 3	feb t-2
%		col 13	dec t-2
%		col 14  jan t-1
%		col 26  jan t
%         col 37  dec t

% E - monthly temperature array

% L1 - logical array specifying rows (years) of D  and E for analysis.
%	L1 is a cv with entries 1 for years to be included in the analysis,
%	and 0 for years not to be included.

% endmo - ending month of tree-ring growth year.  1=jan, 12=dec.
%	So if want to gear analysis to growth year ending in august of
%	year t, set endmo=8.

%********************  CAUTIONS ******************************

% D and E should have been built so that cols 26-37 correspond to jan-dec of
%   year t.  
% D and E should not have missing values over the year range specified by
%   L1

%********************  OUTPUT   *********************************


% Means and standard devs for unwhitened climate series

% Bar graph of corr coef between whitened tri and whitened temp and
% precip for each month (24 bars)

% Bar graph ratio of residual variance to original variance for each
% of 24 clim variables.  Residuals from the multivariate regresssions 
% against all other mothly climate variables.

% Bar graph of AR order, each monthly climate series
% Bar graph of ratio of resid/orig variance from the AR whitening

% The ending month of the monthly grouping varies down the rows.
% The number of months in the grouping varies across the cols.

% Example:  assume endmo=8.

% Col 1 contains r between tree and one-month averages
%	col 1, row 1 is for August of year t
%	col 1, row 2 is for July   of year t
%  col 1,row 12 is for Sept   of year t-1
% Col 2 contains r between tree and two-month averages
%	col 2, row 1 is for July-Aug, year t
%	col 2, row 2 is for June-July, year t
% 	col 2,row 12 is for Aug t-1 thu Sept t-1
% Col 12 contains r between tree and 12-month averages
%	col 12, row  1 is for Sept t-1 thru Aug t
%	col 12, row  2 is for Aug  t1  thru July t
%  col 12, row 12 is for Oct t-2 thru Sept t-1

%*********  PREALLOCATE *********************************************

R = zeros (12);  

S5=zeros(sum(L1),12);
S5T=zeros(sum(L1),12);

lastmo=37- (12-endmo);  % tie ending month to column of D

for k1=1:1;  % loop for each number of months in seasonal grouping
	top=lastmo-(k1-1):lastmo;  % rv.  Say [33] for k1 = 1,endmo=8
	S1 = top(ones(12,1),:); % dupe rows of top
	S2 = (0:11)';
	S3 = S2 (:,ones(k1,1)); % dupe cols of S2
	S4 = S1-S3;  % say, [33 32 31 ...]'

	for k2 = 1:12;  % sum the clim variable over months
		if k1==1; % special case for vector
			S5(:,k2) = D(L1,S4(k2,:));
			S5T(:,k2)=E(L1,S4(k2,:));
		else 
			error('K1 MUST EQUAL 1');
		end ;  % of if loop
	end;  % of for loop k2

	xM = mean([S5 S5T]); % means and variances of monthly ppt, temp
	xS = std([S5 S5T]);
	xM = xM(24:-1:1);
	xS = xS(24:-1:1);

	disp(k1)

end; % of loop for k1
 
