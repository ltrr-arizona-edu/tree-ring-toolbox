function [y1,y2]= meanrat2(B1,jp,mr,md,nomr,mss,P)
%
% Apply mean ratios or mean departures to estimate missing monthly
% PPT or T
%
% D. Meko 10-7-94
%   Last modified 6-23-96 to allow NaN as missing
%
%
%***************   INPUT ARGS  **************
%
% B1 (m1 x n1)r   Data matrix as defined in meanrat1.m (see notes)
% jp (1 x np)i   Sequence number of each predictor series, in order
%	of priority to be considered as the predictor.  
%	Say B1 has 5 cols and jp = [5 2 3 4].  Then col 1 of B1 is the
%	key series and cols 2-5 are data for series 5,2,3,4 as stored
%	in the calling function zest.  This order
%	was set before running the calling function zest, and was
%	stored in the mat file A3 used by zest.m or pest.m.
% mr (1 x np)r   Mean ratios as obtained from meanrat1.m
%	for each of np possible predictors
% md (1 x np)r   Mean departures ....
% nomr (1 x np)L  1=calc of mr or md was not possible for this series
%	0= calc was possible
% mss (1 x 1)r   missing-value code
% P (1 x 1)L    1=PPT,  0=T
%
%
%**************  OUTPUT ARGS  ****************
%
% y1 (m1 x 1)r   B1(:,1) with missing data replaced by estimates
% y2 (m1 x 1)i   sequence number of predictor stations used for the
%	estimates
%
%
%*************   NOTES    ***************
%
% B1 contains PPT or T data for the key station (col 1) and other
% stations specified as possible predictors.  Like meanrt1.m, meanrt2.m
% works on a specific month (say, Jan).  Calling function zest.m loops
% over all 12 months, and has an outer loop over all key stations.
%
% Checks that estimation is needed and possible.  Estimation not needed
% if no missing data at key station this month.  Then y1 will be set
% to the original data, and y2 will be set to vector of zeros.
%
% Estimation might not be possible for several reasons.
%  1. calc of mean ratio/departure might not have been possible for
%     any of the specified potential predictors.  Possible cause
%     is no overlap of nonmissing data;  or, if PPT, only zeros in
%     the overlap.  These cases are covered by nomr(n)
%  2. For the year in question, nonmissing data might not exist at
%     any of possible predictor stations.
%
% Estimate is set to zero when all of these conditions are true:
% (1) at least one predictor station has zero ppt, and rest have zero or
%     missing values
% (2) mean-ration estimation was otherwise not possible, for example, if
%     the predictor station with zero ppt had no overlap of data with the
%     key station in other years with which to calc the mean ration, or
%     if data in the overlap was zero, so that an infinite mean ration 
%     would result.

[m1,n1]=size(B1);
z = B1(:,1);  % key station data
y1 = z;  % initialize output vector as unmodified key data
y2 = zeros(length(y1),1);  % initialize predictor ident matrix to 0
V = B1(:,2:n1);  % predictor-station data

% Check that there is a need to estimate
L1p = z ~= mss & ~isnan(z);   % data present at key station (cv)
L1m = ~L1p;      % data missing at key station (cv)

L2 = ~isnan(V) & V ~= mss;  % nonmissing elements of predictor matrix

% Expand rv nomr into a matrix same size as V
NOMR = nomr(ones(m1,1),:);

% Point to elements nonmissing and estimatable
L3 =  L2 & ~NOMR;

% Point to years where any data with above conditions
L4 = (any(L3'))';   %  cv, years with data nonmissing and estimatable

L5 = L1m & L4;  % key data missing, predictor data present and estmbl


% The estimation loop
if  sum(L1m) == 0   |  sum(L5) == 0 ; % no need to estimate or not poss
   % no action needed, except possibly special case at end of code
else
   % Pull subset of years estm to be done for, and cols of 
   % valid predictors (those with nomr(n)=0)
   z1 = z(L5);  % key series
   mz1 = length(z1);  % this many values will be estimated
   V1 = V(L5,~nomr);

   % Pull segment of mr, md vectors for valid predictors
   mr1 = mr(~nomr);
   md1 = md(~nomr);
    
   % Get predictor sequence numbers of the valid predictors
   jp1 = jp(~nomr);
   nnp = length(jp1);  % This many predictors in pool

   % Logical operations to make logical matrix 
   % whose columns will tell wich years in z1 will be estimated
   % using first predictor series, second series, etc
   LL1 = (L2(:,~nomr))';  % data present in predictor elements;  transpose so that
     % cols corresp to years
   [mLL1,nLL1]=size(LL1);

   if mLL1==1 ;  % handle special case of a 1-row matrix
   	LL2 = LL1;
   else
	LL2 = cumsum(LL1);
   end

   LL2 = LL2 == 1;
   LL3 = LL2 & LL1;
   LL4 = LL3';   
   LL4 = LL4(L5,:);
   % Col 1 of LL4 has 1s in years that jp1(1) is to be predictor
   % Col 2 of LL4 ...                  jp1(2)
   % etc

   % Initialize data vector with unmodified values, and begin estimation
   yy1 = z1;
   yy2 = (zeros(length(z1),1));
   for n = 1:nnp;  % loop for each valid possible predictor
	jp1n = jp1(n);  % seq number of this predictor
	LLL = LL4(:,n);
	nL = length (LLL);
	sL = sum(LLL);
	if any (LLL);  % At least one year can be estim from this predictor
	   x=V1(LLL,n);  % cv of predictor values for years to be used
	   if P;  % PPT, not T
	      est = x * mr1(n);
	   else;  % T
	      est = x + md1(n);
	   end
	   yy1(LLL) = est;
	   yy2(LLL) = jp1n(ones(sL,1),:);
	else
	end
   end
 
   y1(L5) = yy1;
   y2(L5) = yy2;

end

if P==1
   % Special case of P data with
   %  1. key data missing
   %  2. only non-missing predictor data zero
   %  3. other screening determined not possible to estimate, say because
   %     no acceptable data for computing mean ratios.  
   min1 = -1;
   L6 = V == 0;  % predictor data zero
   L7 = ~L2 | L6 ;  % predictor data either missing or zero
   % Mark all years not estimated by above method, but
   % where the only nonmissing predictor data is zero
   L8 =  L1m &  ~L5   &    (all(L7'))'  &  (any(L6'))';
   nL8 = sum(L8);   
   y1(L8)=zeros (nL8,1);
   y2(L8)=min1(ones(nL8,1),:);
else
end
	
