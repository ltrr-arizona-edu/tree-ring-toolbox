function xe = dayest2(xk,xp,xmiss,xeold)
% Estimate missing daily ppt for each of 366 julian days for a station
%
% D. Meko 11-12-93
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
% %***************** USE **************
%
% Called by dayfill.m , which fills in a matrix of
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
% Rules for various cases:
%
%	1. xeold missing, xp zero ----  substiture zero
%	2. xeold missing, xp present and non-zero -- median ratio estimation
%	3. xeold missing, xp missing ---- no action
%
%	Each julian day treated separately.  Only observations with nonzero
%  data for both xk and xp are used in the subsample for calc of median ratio.  
%  A ratio series is computed of  xk/xp, using as a sample all acceptable
%  observations in the 7-day period centered on the julian day.  The median
%  of this ratio series is then multiplied times the values in the 
%  predictor series (xp) to get the estimates for xe.
%
%
%
%**************** USER-WRITTEN FUNCTIONS NEEDED 
%
% medday.m -- estimation by median ratio method
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
	error('length of xk should be a multiple of 366 ')
end
nyr = mx/366; % number of years of data

if mx ~= [length(xeold)  length(xp)];
	error('xp, xeold, xk must be same length');
end

L=zeros(mx,4); % logical vector
JT1=zeros(366,1);
JT2=zeros(366,nyr);
LS=zeros(nyr*366,1); 

xe = xeold;  % initialize xe

% Fill cols of L with pointers to rows of xk, xp, xeold satisfying
% specified conditions 
L(:,1) =  xk~=xmiss & xp ~= xmiss; % xp and xk not missing 
L(:,2) = xe == xmiss & xp~=xmiss ; 
L(:,3)= xk~=0 & xp~=0; % xk and xp non-zero
	% subsample to be estimated
L(:,4) = xe == xmiss & xp==0;  % sample for substituting zero

% Don't go further if no missing values in xeold with corresponding nonmissing
% in xp
if ~any(L(:,2))
	disp('No work to do.');
	disp('xeold has no missing values, or ');
	disp('none that go with nonmissing in xp');
	return
end

% Make vector of Julian days, stacked nyr times
jd1 = (1:366)';
jd1 = jd1(:,ones(nyr,1));  % dupe vector to nyr cols
jd = reshape (jd1,366*nyr,1);  % reshape to col vector


J=zeros(366,1);
J(1:4)=ones(4,1);
J(364:366)=ones(3,1);
JT=toeplitz(J);

% Before entering day loop, replace all missing values with zero if
% corresp value for predictor series is zero
if any(L(:,4));
	xe(L(:,4))=zeros(sum(L(:,4)),1);
end

% Loop through each julian day.  Samples for median ratios are
% from 7-day periods centered on julian day
for i=1:366 ;
	JT1=JT(:,i);
	JT2=JT1(:,ones(nyr,1));
	LS=reshape(JT2,nyr*366,1);   % points to 7-day subsets, this julian day

	% Pointer to rows in xe that need to be estimated by median ratio
	L2 = jd==i; %  this julian day
	L3 = L2 & L(:,2) & xp~=0; %  this jd, missing in xe, present in xp,
		% and not zero in predictor series
	L4 = L(:,1) & L(:,3) & LS;  % subsample for med ratio calcs

	if any(L3);  %  At least one missing value needs to be estimated
		if any(L4);  % At least one case with nonzero nonmissing data for
			% both xk and xp.  Otherwise, can't form ratio series.
			xe(L3)=medday(xk(L4),xp(L4),xp(L3)); % estimate
		end
	end
	disp(['Finished Julian-day ',int2str(i)]);

end
