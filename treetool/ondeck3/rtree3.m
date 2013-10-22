function r = rtree3(X,nms,yrs)
% rtree3: correlation coefficient between tree indices
% r = rtree3(X,nms,yrs);
% Last Revised 9-2-99
%
% Correlation coefficients between tree indices, as needed by stats2.m and treei
%
% INPUTS :  X (? x 1) -  Tree indices
% 	    nms (? x 1) -  Names of the tree indices
% 	    yrs (? x 3 ) -  Year and index info on X
%	   
% OUTPUTS :  r (? x 5) -  Correlations and related info 
%		Column 1 - Correlation coeff
%		Column 2 - Sample size.
%		Column 3 - Tree sequence number of first tree
%		Coumns 4 - Tree sequence number of second tree
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%*** NOTES
% Any masking out of trees from the correlations is assumed
% to have been done before call to this function. In other words, 
% correlations are generated for all series in the string matrix nms.
%
%____________________________________________________________


[ns,dum]=size(nms); % ns is number of trees
r = zeros(ns*(ns-1)/2,4);  % initialize

if(ns==1); % just one tree, return r empty
	r=[];
	return
end

row = 0;  %initialize row of r

for  k1 = 1:(ns-1); % loop over key trees
	% Get start and end row index in X of key series 
	i1 = yrs(k1,3);
	i2 = i1 + (yrs(k1,2) -yrs(k1,1));
	x = X(i1:i2); % get key series, full length
	yrx = (yrs(k1,1):yrs(k1,2))'; % year vector for key series
	
	for k2 = (k1+1):ns; % loop over comparison series
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
				r(row,3)=k1;
				r(row,4)=k2;
			else
				r(row,:)=[NaN NaN NaN NaN];
			end
		else
				r(row,:)=[NaN NaN NaN NaN];
		end; % overlap>=20 if
	end; % for k2
end; % for k1

