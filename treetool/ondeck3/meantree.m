function [IT,ITn,Tnms,Tyrs] = meantree(IX,nms,yrs,msk)
%
% DESCRIPTION : [IT,ITn,Tnms,Tyrs] = meantree(IX,nms,yrs,msk)
% Converts standard core indices, as culled by the mask from 
% BLDMASK.M, into tree indices by averaging all available core 
% indices for each tree in each year.
%
% INPUTS  :  IX (? x 1)   - Strung out vector of core indices
%	     nms (ns x 6) - Names for the cores
%	     yrs (ns x 3) - Years and starting row subscripts 
%			    for cores within IX
%       msk (? x 1)L  mask: 1=use this core, 0 = mask out this core
%
% OUTPUTS :  IT (? x 1)   - Strung out vector of tree indices
% 	     ITn (? x 1)  - Contain averaging info. Same length
%			    as IT.
%	     Tnms (? x 6) - Tree names
%	     Tyrs (? x 3) - Years and starting row indices for 
%			    trees in IT 
%______________________________________________________________
%
% UW FUNCTIONS CALLED
% treenum
% maskind

% Get parameters from treenum.M.  
[I,Tnms,n]=treenum(nms,msk);
  % n is number of tree-indices that will be produced
  % I holds the sequential core numbers (col 1) and the 
  %   sequential tree numbers (col 2).  Entries in col 2 will
  %   go from 1 to n.
  % T holds the 6-character (blank-filled) tree ids for the n trees


% Figure out the required length of the strung-out column vectors
% IT and ITn, and initialize as NaN
% - Loop over trees, finding earliest start year of any core and
%   latest end year of any core, not counting leading or trailing
%   NaNs.  Will need to first build a dummy ? x 10 NaN matrix,
%   initializing for each tree. Then pull the data from IX and
%   store in the dummy matrix.  Then check for earliest and latest
%   year in this matrix with not all values in row NaN
a = NaN;
L1 = I(:,2)~=0;  % points to active rows of yrs
YRS = yrs;
yrs1 = yrs(L1,:); % subset of rows of yrs
[m1,n1]=size(yrs1);  % m1 is number of active cores these trees
for j = 1:m1; % loop over cores
	years = yrs1(j,:);
	p1 = years(3); % start row of this core index in IX
	p2 = p1+ years(2)-years(1); % end row
	x = IX(p1:p2);  % pull the index
	% find first and last valid (not NaN) value
	L2 = ~isnan(x);
	f1 = find(L2);
	r = [min(f1)  max(f1)]; % rows of first and last valid data in x
   % Make an adjusted index into IX
  
	y1 = years(1)+r(1)-1;  % first year of valid data
	y2 = years(2)-(length(x)-r(2));  % last year of valid data
	p3 = p1 + r(1) - 1;   % row index in IX of year y1
	p4 = p2 - (length(x)-r(2)); %row index in IX of year y2
	yrs1(j,:)=[y1 y2 p3];
end
YRS(L1,:)=yrs1;  % like yrs, but with revised years and row index
	% to account for lopped off NaNs on ends
	

% Size and allocate strung-out vectors IT and ITn.  Initialize as
% col vectors filled with NaN
nsize = 0; % initialize length  of IT and ITn to zero
for i = 1:n;  % loop over trees
	L3 = I(:,2)==i; % logical pointer to rows in YRS with years info
		% for this tree's cores
	nc(i) = sum(L3);  % number of cores for this tree
	A1 = YRS(L3,:);  % years row vector info for this tree's cores
	A2 = [min(A1(:,1))  max(A1(:,2))]; % start, end year for this tree
	Tyrs(i,1:2)=A2;
	nsize = nsize+A2(2)-A2(1)+1;  % extend needed length for IT and ITn
end


% Compute row index of first tree-index value for each tree in IT
D1 = Tyrs(:,2)-Tyrs(:,1)+1; % number of years for each tree index
D2 = cumsum(D1);
D2(length(D2))=[];
D2=D2+1;
D2=[1;D2];
Tyrs(:,3)=D2;

% Initialize the strung-out column vectors as NaNs
IT = a(ones(nsize,1),1);
ITn = IT;



%****************  COMPUTATION OF TREE INDICES *******
% 
% Now have key information:
%   IX ,yrs:  core-indices vector and years pointer
%   YRS -- pointer to core-index vector to avoid the leading
%		or trailing NaNs, which will be there if the IX is
%		residual index, or if have truncated parts of ring-width
%		series in curve fitting
%   IT, ITn:  the col vector (initialized to NaNs) in which
%		the tree indices will be stored
%   Tyrs: the years and starting-row pointer to IT and ITn telling
%		which slots the individual tree indices will be stored




for i=1:n; % Loop over trees

	% Size the temporary ? x 10 data-storage matrix and
	% Initialize the matrix as NaN
	U=Tyrs(i,:);
	nrows = U(2)-U(1)+1;
	Z = a(ones(nrows,1),ones(10,1));

	% Pull the row-index info into IX for the cores for this tree
	L3 = I(:,2)==i; % logical pointer to rows in YRS with years info
	G = YRS(L3,:);
	[mG,nG]=size(G);  % mG is number of cores for this tree
	I4 = G(:,1) - U(1) + 1; % starting slots in Z for core indices
	I5 = I4+G(:,2)-G(:,1); % ending slots

	% Loop over the cores to be used for this tree
	for j = 1:mG;
		i4 = I4(j);
		i5=I5(j);
		k5 = G(j,3); % starting row index of the core data in IX
		k6 = k5 + (G(j,2)-G(j,1)); % ending row index
		z = IX(k5:k6);
		Z(i4:i5,j)= z;
	end

	% Compute vector of sample-size (number of cores in each year)
	L5 = ~isnan(Z);
	tn = (sum(L5'))';


	% Substitute zeros for NaNs in the sub-matrix
	L6 = ~L5;
	numz = sum(sum(L6));
	zee = zeros(numz,1);
	Z(L6)=zee;

	% Compute the tree mean, a col vector


	t = (sum(Z'))' ./ tn;


	% Put tree mean and its sample size in slots of IT and ITn
	k7=U(3); % starting index
	k8 = k7 + U(2)-U(1); % ending index
	IT(k7:k8) = t;
	ITn(k7:k8) = tn;
end


% Put NaNs in IT where sample size is zero (would replace zeros)
L8 = ITn == 0;
zee = sum(L8);
if zee>0;  IT(L8) = a(ones(zee,1),:); end;
