function [mr,md,nomr,ss,rr]=meanrat1(B1,np,mss,P)
%
% mean ratios to be used in zest.m and meanrat2.m
%
%
%*************  INPUT ARGS **********************
%
% B1 (m1 x np+1)r   matrix of monthly PPT or T;  each row a year
%	col 1 is key (predictand) series
%	cols 2, 3,... are predictor series 
% np (1 x 1)i   number of stations in pool of possible predictors
% mss (1 x 1)i  numeric missing-value code [-999]
% P (1 x 1)L   1 if PPT, 0 if T
%
%**************   OUTPUR ARGS  *****************8
%
% mr (1 x np)r   mean ratios, ratio of mean PPT at key station to 
%		mean PPT at each of np possible predictor stations;
%		ratios based on data for same period at both key and
%		predictor station.  Filled with mss if nomr(i)=1.
% md (1 x np)r   mean departure, difference of mean T at key and 
%		predictor stations.
% nomr (1 x np)l   1=calc of mean-ratio or median ratio not possible
% ss (1 x np)i   sample size for mean ratio or mean departure
% rr( 1 x np)r   correl coeff between key series and each predictor, 
%	for period mean ratios computed on;  NaN returned if sample
%	size fewer than 5
%
%
%************* NOTES ***********************
%
% Works on a block of data for a particular month (say, Jan).  Typically,
% the block covers all years of data at the key or predictand series.
% Computes both the mean ratios and mean departures, regardless if data
% is PPT or T.  Function meanrt2.m will use the correct version.
%
% Mean ratio or mean departure cannot be computed  

z = B1(:,1);   % data at key station
V = B1(:,2:np+1);  % data at predictor stations


ss = zeros(1,np);  % allocate for sample size for mean ratio
notanum = nan;
rr = notanum(:,ones(np,1));

for n = 1:np;  % Loop over predictors
	v = V(:,n);  % cv of data for this predictor
	% Make pointer to years that will be used to calc mean ratios
	% or departures
	L1 =  (~isnan(z) & z ~= mss) & (~isnan(v) & v ~= mss) ; 

	if sum(L1) == 0;  % no year with data at both key and predictor
	   nomr(n)=1;
	else
	   nomr(n) = 0;
	end


	if ~nomr(n) ; % if at least one year with data at key and predic
		v1 = v(L1);  % subset of predictors for years with  for those years
		z1 = z(L1);
		if P & ((max(v1) == 0) |  (max(z)==0)); % neither key nor predict
			% with any non-zero precip data for those years with nonmissing
			% data at key and predictor
			nomr(n) = 1;
			ss(n)=nan;
		end
	end


	if ~nomr(n);  % Ok, mean-ratio or mean departures can be 
		% calculated for this series.
		mr(n) = mean(z1) / mean (v1);
		md(n) = mean(z1) - mean (v1);
		ss(n)= length(z1);
		if length(z1)>=5;
			rrr = corrcoef(z1,v1);
			rr(n) = rrr(1,2);
		end
	else;
		mr(n) = nan;
		md(n) = nan;
		ss(n) = nan;
		rr(n) = nan;
	end
		
end
	

