function xe = dayest(xk,xp,xmiss,xeold)
% Estimate missing daily ppt for each of 366 julian days
%
% D. Meko 11-8-93
%
%
% xk (mx x 1) daily ppt for 366 julian days in each of unk # of yrs
%		at key station -- the station for which data will be estimated
% xp (mx x 1) daily ppt likewise for a 'predictor' station, the 
%		station used to estimate the missing values
% xmiss (1 x 1) missing-value code -- numeric
% xeold (mx x 1) vector similar to xk, but which may already have
%		had some missing values estimated -- for example, in a 
%		previous run of this function using a different predictor stn
%
%
% xe (mx x 1) updated daily ppt at the key station; includes the
%		newly estimated values using xp
%
%
%***************** USE **************
%
% Function under dayfill.m.  One of set of functions for filling in a
% matrix of daily ppt values at several stations by estimating missing
% data.  See dayfill.m for method.
%
% dayest.m handles estimation for a specified pair of series -- one
% the 'key' station, the other the 'predictor' station.  A loop 
% goes through all 366 Julian days.  
%
%
%**************** NOTES **************
%
% Procedure for various cases:
%
%	1. xk missing, xp zero ----  xe=x=0
%	2. xk missing, xp present and non-zero
%		2.1  No non-missing subsample of xk --- xe=xp
%		2.2  Subsamples present ----  xe estimated by nonexceed prob
%	3. xk missing, xp missing ---- no action
%
%
%
%**************** USER-WRITTEN FUNCTIONS NEEDED 
%
% freqest1.m -- estimation using nonexceedance probabilities
%
%
%*************** GLOBALS ***************
%
%  none
%
%
%*************  PREALLOCATE AND BEGIN *********


mx=length(xk);  % if nyr years, will be 366 x nyr
if (rem(mx,366)~=0),
	error('xk should have a multiple of 366 years')
end
nyr = mx/366; % number of years of data

if mx ~= [length(xeold)  length(xp)];
	error('xp, xeold, xk must be same length');
end

L=zeros(mx,5); % logical vector
xe = xeold;  % initialize xe

% Make vector of julian days, stacked nyr times
jd1 = (1:366)';
jd1 = jd1(:,ones(nyr,1));  % dupe vector to nyr cols
jd = reshape (jd1,366*nyr,1);  % reshape to col vector

% Fill cols of L with pointers to rows of xk, xp, xeold satisfying
% specified conditions
L(:,1) = xp ~= xmiss & xp>0;  % xp not missing or zero
L(:,2) = xk ~= xmiss & xk >0;  % xk not missing or zero
L(:,3) = L(:,1) & L(:,2); % xp and xk not missing or zero
	% subsample for cum freq distributions
L(:,4) = xe == xmiss & (xp~=xmiss & xp~=0);
	% subsample that will use exceedance probs for estimation
L(:,5) = xe == xmiss & xp == 0; 
	% subsample whose missing values will be replaced with zero




% Loop through each julian day.  Samples for cum freq distribs will
% be three-day periods centered on the Julian day.
for i=1:366;
	% Pull 5-day samples
	if i==1; % Special case for Julian day 1
		L3=jd==365 | jd==366 | jd==1 | jd==2 | jd==3;
	elseif i==2;
		L3= jd==366 | jd==1 | jd==2 | jd==3 | jd==4;
	elseif i==365;
		L3= jd==363 | jd==364 | jd==365 | jd==366 | jd==1;
	elseif i==366;
		L3=jd==364 | jd==365 | jd==366 | jd==1 | jd==2;
	else
		L3= jd== (i-2) | jd==(i-1) | jd==i | jd==(i+1) | jd== (i+2);
	end;  % of if for special cases days 1 and 366
	L3p = L3 & L(:,1); % points to years with xp not missing or zero
	L3k = L3 & L(:,2);  % same for xk
	L3pk= L3 & L(:,3); % points to years in which both xk and xp
		% are non-missing and non-zero.

	iopt=0;
	if ~any(L3p) | ~any(L3k);
		disp('No years with non-missing and non-zero at one stn');
		pause(3);
		iopt=1;
		% Cannot use exceed prob method in this case.  So merely 
		% substitute the value for the predictor series as the 
		% value for the key series
	end

	% Pointer to rows in xe that need to be estimated using exceedance
	% probs.
	L4 = jd==i; % point to this julian day
	L5 = L4 & L(:,4); %  this jd, missing in xe, present and nonzero
		% in xp

	% Special case of direct substitution
	if iopt==1;
		xe(L5)=xp(L5);
	end

	if (any(L5) & iopt==0);  % Go thru exceedance-prob estimation only if need to
		xe(L5)=freqest1(xp(L3p),xk(L3k),xp(L5)); % estimate
	 end;  % of if any(L5)

	% If xe is missing and xp is zero, set xe to zero
	L6 = L4 & L(:,5);
	if any(L6);
		xe(L6) = zeros(sum(L6),1);
	end
	disp(['Finished Julian-day ',int2str(i)]);
end


