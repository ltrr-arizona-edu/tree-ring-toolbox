function n=treecnt1(nms,yrs,t,yr1);
%n=treecnt1(nms,yrs,t,<yr1>);
%
% Number of trees at a site in specified years
%
% Meko 3-1-97
%
%
%****************** IN  *********************************
%
% nms (:,8)s  site-tree-core ids, in sov storage format
% yrs (:,3)i  start year, end year, row index of start 
%		year for each core.  row index not used in this function.
% t (:,1)i  vector of target years 
% yr1 (1 x 1)i  optional, specifies that not only must
%		the tree exist in years specified by t, but also
%		must still be there in year yr1.  yr1 might be a 
% 		recent year (e.g., 1968), so that you want to know
%		if the tree's data could be calibrated directly 
% 		with climate data
%
%********************* OUT **************************
%
% n (:,1)i  number of trees in years specified by t
%
%
%***************   FUNCTIONS CALLED ******************
%
% treenum.m  analyzes string matrix of names to find unique trees
%
%*********************** NOTES ***********************
%
% See notes in treenum.m for how the ids are examined to arrive
% at a tree count.
%


[nyears,nt]=size(t);	
if nt~=1;
	error('t must be cv')
end

a=NaN;
n=a(ones(nyears,1),:); 

% Loop over years in t
for  i = 1:nyears;
	year = t(i); % target year
	
	% Is the target year in t enclosed by time period of cores?
	L1=yrs(:,1)<= year & yrs(:,2)>=year;

	% Recent year constraint?
	if nargin==4;
		L2=yrs(:,2)>=yr1;
		L1=L1 & L2;
	end

	if sum(L1)>1;
		nms1=nms(L1,:)  % pointer to selected cores

		% Form a moot core mask needed as an input arg to treenum.m
		m1=size(nms1,1);
		cmask=ones(m1,1);

		% Calculate the number of trees represented in the selected cores
		[I,T,ntrees]=treenum(nms1,cmask);  % I, T not used here
	elseif sum(L1)==1;
		ntrees=1;
	else
		ntrees=0;
	end
	 % Store result
	n(i)=ntrees;
end
