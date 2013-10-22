function r = rtree1a(X,nms,yrs,mask)

% DESCRIPTION : r = rtree1(X,nms,yrs)
% This function finds the between-core correlation coefficients 
% for pairs of cores and summarizes results according to whether
% cores from same or different tree
%
% INPUTS :  X (? x 1) -  Core indices
% 	    nms (? x 1) -  Names of the core indices
% 	    yrs (? x 3 ) -  Year and index info on IX
%      mask (? x 1) -- masks out correlations for summary of between
%			and within-tree correlations
%
% OUTPUTS :  r (? x 5) -  Correlation matrix. 
%		Column 1 - Correlation coeffs.
%		Column 2 - Sample size.
%		Column 3 - 1 if same tree, 2 if not.
%		Coumns 4 & 5 - Sequence numbers of Core 
%			ID's used to find the corr. coef; corresponds to
%			row number of the core in original nms.
%
%____________________________________________________________


[mm1,nn1]=size(mask);
cdum = ones(mm1,1);
[I,Tdum,ndum]=treenum(nms,cdum);
t = I(:,2); % sequential tree number for each core

[ns,dum]=size(nms); % ns is number of trees
r = zeros(ns*(ns-1)/2,5);  % initialize

if(ns==1); % just one tree, return r empty
	r=[];
	return
end

row = 0;  %initialize row of r


for  k1 = 1:(ns-1); % loop over cores
	% Get start and end row index in X of key series 
	i1 = yrs(k1,3);
	i2 = i1 + (yrs(k1,2) -yrs(k1,1));
	x = X(i1:i2); % get key series, full length
	yrx = (yrs(k1,1):yrs(k1,2))'; % year vector for key series
	jtree=t(k1);
	
	for k2 = (k1+1):ns; % loop over comparison series
		ktree=t(k2);
		if jtree==ktree,
			kk = 1;
		else
			kk=2;
		end
		row=row+1; % increment row counter for r
		% find overlap
		yr1 = max(yrs(k1,1),yrs(k2,1)); % latest start year
		yr2 = min(yrs(k1,2),yrs(k2,2)); % earliest end year
		overlap = yr2-yr1+1;
		% Require minimum of 20 years of overlap
		if overlap >= 20,
			% Get start,end row index in X for comparison series
			i3 = yrs(k2,3);
			i4 = i3+(yrs(k2,2)-yrs(k2,1));
			y = X (i3:i4); % get comparison series, full length
			yry = (yrs(k2,1):yrs(k2,2))'; % year vector for comp.series

			% Make logical pointers to rows of x,y in overlap period, then
			% pull out the common-period data
			Lx1 = yrx >= yr1 & yrx <= yr2;
			Ly1 = yry >= yr1 & yry <= yr2;
			x1=x(Lx1);
			y1 = y(Ly1);

			% Remove rows with any NaNs
			z = [x1 y1];
			L3 = isnan(z);
			L3sum = (sum(L3'))';
			L4 = L3sum>0;
			x1(L4)=[];
			y1(L4)=[];
			n1=length(x1);  % number of years correl based on
			if n1>=20;  % another check on min sample size
				dum = corrcoef([x1 y1]);
				r(row,1)=dum(1,2);
				r(row,2)=n1;
				r(row,3)=kk;
				r(row,4)=k1;
				r(row,5)=k2;
			else
				r(row,:)=[NaN NaN NaN NaN NaN];
			end
		else
				r(row,:)=[NaN NaN NaN NaN NaN];
		end; % overlap>=20 if
	end; % for k2
end; % for k1





