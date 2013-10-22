function R = r4sov(X,YRS,mlap,p,kopt)
% R = r4sov(X,YRS,mlap,p,kopt)
%
% Correlations between pairs of time series stored in a 
% sov.
%
% Meko 2-19-97
%
%************* IN ****************************************** 
% X (mX x 1)r strung out vector of time series, one after another
% YRS (nums x 3)i start year, end year, and row index of starting 
%		value of each time series in X.  nums equals the number
%		of time series  stored in X
% mlap (1 x 1)i minimum overlap to compute a correl coef on. 
%		Correl returned as NaN if not satisfied
% p (1 x 2)i  start and end year of period for which correlations
%		desired.  For example, might have a set of tree ring indices,
%		and want correlations between cores for some recent calibration
%		period. Period of analysis will be restricted to period p
%		or whatever subset of years of p with valid data available. See
%		kopt(1)
% kopt (1 x 1)i options
%		kopt(1)	==1 no restriction on when overlap must occur
%				==2 overlap must occur in period specified by p
%
%********************** OUT **************************************
%
% R (? x 5) -  Correlation matrix. 
%		Column 1 - Correlation coeffs.
%		Column 2 - Sample size.
%		Coumns 3 & 4 - Sequence numbers of time series telling
%			which pair of series the correlation coef applies toCore 
%			ID's used to find the corr. coef; corresponds to
%			row number of the core in original nms.
%
%____________________________________________________________

[nums,ntemp]=size(YRS); %nums is number of time series
if ntemp~=3,
	error('YRS col size must be 3');
end

if(nums==1); % just one time series, return R empty
	R=[];
	return
end

R = zeros(nums*(nums-1)/2,4);  % initialize

row = 0;  %initialize row of R

for  k1 = 1:(nums-1); % loop over time series to be first of pair (key)
	% Get start and end row index in X of key series 
	kstop=0;
	i1 = YRS(k1,3);
	i2 = i1 + (YRS(k1,2) -YRS(k1,1));
	x = X(i1:i2); % get key series, full length
	yrx = (YRS(k1,1):YRS(k1,2))'; % year vector for key series
	
	% optionally compute number of years of valid data in key
	% series in period p
	if kopt(1)==2;
		L5=yrx>=p(1) & yrx<=p(2);
		if sum(L5)<mlap,
			kstop=1;
		else; % grab the segment of x in period p
			yrx=yrx(L5);
			x=x(L5);
			L6=~isnan(x);
			if sum(L6)<mlap
				kstop=1;
			else
				x=x(L6);
				yrx=yrx(L6);
			end
			d1=diff(yrx);
			if ~all(d1==1);
				error(['Noncontinuous years, key series ',int2str(k1)]);
			end
		end
	end; % of if kopt(1)==2

	if kopt(1)==2 & kstop==1; % not enough valid overlap, set 
	% corrs to NaN
			% no action needed;

	else
		
		for k2 = (k1+1):nums; % loop over comparison series
			row=row+1; % increment row counter for R
			% find overlap
			yr1 = max(YRS(k1,1),YRS(k2,1)); % latest start year
			yr2 = min(YRS(k1,2),YRS(k2,2)); % earliest end year
			overlap = yr2-yr1+1;
			% Require minimum of 20 years of overlap
			if overlap >= mlap,
				% Get start,end row index in X for comparison series
				i3 = YRS(k2,3);
				i4 = i3+(YRS(k2,2)-YRS(k2,1));
				y = X (i3:i4); % get comparison series, full length
				yry = (YRS(k2,1):YRS(k2,2))'; % year vector for comp.series

				% Make logical pointers to rows of x,y in overlap period, then
				% pull out the common-period data
				Lx1 = yrx >= yr1 & yrx <= yr2;
				Ly1 = yry >= yr1 & yry <= yr2;
				x1=x(Lx1);
				y1 = y(Ly1);

				% Remove rows with NaN in either series
				z = [x1 y1];
				L3 = isnan(z);
				L3sum = (sum(L3'))';
				L4 = L3sum>0;
				x1(L4)=[];
				y1(L4)=[];
				n1=length(x1);  % number of years correl based on
				if n1>=mlap;  % another check on min sample size
					dum = corrcoef([x1 y1]);
					R(row,1)=dum(1,2);
					R(row,2)=n1;
					 R(row,3)=k1;
					R(row,4)=k2;
				else
					R(row,:)=[NaN NaN NaN NaN ];
				end
			else
					R(row,:)=[NaN NaN NaN NaN ];
			end; % overlap>=20 if
		end; % for k2
	end; % of if kopt(1)==2 & kstop==1
end; % for k1

% Get rid of all-NaN rows
L5=isnan(R);
L6=(all(L5'))';
R(L6,:)=[];
