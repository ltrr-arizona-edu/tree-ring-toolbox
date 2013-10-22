%  medcorr.m      % median correlation 

% D. Meko , Sept 25, 1991

% Computes median and maximum correlation of tree-ring series with
% specified climate variables from array B1.

% First used along with near1.m to compute correlations between 
% each tree-ring chron and all hcn ppt stations within a 
% 110 km search radius.
% Near.m provides the distance computations, and 
% gives the subscripts of the approp clim variables from B1.

% All tree-ring series (cols of A1) need not be used in the analysis.
% If search distance in near.m is small, some tree sites may have no
% nearby climate stations.  You will want to run medcorr.m only on the
% subset of tree-ring sites with valid climate-station neighbors.  The
% index array I2 is used to specify which tree-ring sites are to be 
% used.  I2 is an index into cols of A1.  If have valid neighbors for
% all tree sites, I2 is identically the cv (1:n1)', where n1 is number
% of cols in A1.

% Depending on the season of the climate variable, the climate
% series may have to be shifted by one year relative to the tree-rin
% series to offset the most physically reasonable monthly groupings.
% For example, Fall, year n  includes Sept-Nov year n, which most
% logically would be correlated with tree rings for year n+1.  So, if
% the climate series goes from 1906-87 and the trees from 1755-1979,
% the period for analysis should be 1907-79, and the climate series 
% should be shifted back 1 year before the correlation analysis.  A screen
% prompt will handle this shift.

%**********   BEFORE ANYTHING  *************************************

% Run near.m on the files of x,y coordinates to generate arrays D1, D,
% and nD.

% Form the index array I2 by I2=find(nD>0).

% Save D, D1, nD and I2 in dist.mat.   The arrays hold:

% D1 (m1 x ndcols) [I]  The subscripts of the climate stations representing
%	the nearest stations to each of n1 tree-ring sites.  Although D1
% 	has ndcols columns, there may be fewer than ndcols valid subscripts.  For 
%	example, if the search distance option yields only 5 stations, 15
% 	remaining values will be invalid.  These invalid elements are set
% 	to -2.0.

%	The subscripts make sense only relative to existing tree-ring and 
% 	climate arrays.  These arrays were required by near.m

% D	(m1 x ndcols) [R]  The corresponding array of distances between m1
%	tree-ring sites and each climate station.  Invalid elements filled
%	with 2's.

% nD	(m1 x 1) [I]  The number of valid nearest neighbors for each tree
%	ring site.  Less than or equal to ndcols.

% I2 index cv into cols of A1, specifying which of the tree sites have
%	 neighbors within search distance.


%*****************   PRELOADS   ***********************************

%  1. 	dist.mat   (see above)
%  2.		A1 (m1 x n1) [R] the tree-ring index array
%  3.		B1 (m2 x n2) [R] the climate array
%  4.		E1 (3 x 2) [I]  first and last years of 
%			row 1 -  the tree-ring array
%			row 2 -  the climate array
%			row 3 -  the analysis period.



%*****************  INPUT   *******************************************

%  A1, B1, E1, D, D1, nD, I2  -- as described above


%******************   OUTPUT   ***************************************


%  R2 (n1 x 20) [R]  corrln of each tree ring series (row) with up to
%		ndcols climate stations.  Invalid data filled with zeros.
%	   Order matches stn subcripts and distances from D,D1.  Only nD(i)
%		correlations are valid for ith tree-ring site.

%		Can later plot dependence of correlation on distance from tree site
%		by plotting the first nD(i) members of R2 against the first nD(i)
%		members of D, for row i of R2 and D.

%  R4 (n1 x 5) [R] if saved ascii, surfer-readable for plotting a map of
%     the correlation coefs.  Col-1,2 the x,y coords, in km
%	col-3  the number of stations in the search distance
%  col-4  the median correlation
%  col-5  the maximum correlation

%  nD (n1 x 1) [I] how many valid correls in each row of R2.  This vector
%	is actually known previously from near.m.  Near.m points out the
%	stations to be considered in the correlation analysis for each tree
%	ring site.  


%*****************    OTHER VARIABLES  ******************************

% IA years of input tree rings
% IB years of input climate data
% IX years of intended analysis period.
% m1,n1  size of A1
% m2,n2  size of B1
% n4  equiv to length(I2), the number of tree sites valid for analysis
% m3 number of years in intended analysis period.
% R (n1+n2,n1+n2)  correlation matrix of tree rings and climate
% R1 (n1 x 4) [R] summary correln info for each site:
%  	col 1 - the site sequence no.
%  	col 2 - the number of stations (valid in search dist) 
%	col 3 - the median corrln coef of the valid entries
%	col 4 - the maximum corrln
% tree (n1 x 2) x,y coords for each tree-ring site

%*******************   TRIM ARRAYS  **********************************

[ndrows,ndcols]=size(D1);

skip= input('QUICKY DEBUG MODE--NO REFORMING OF A1,B1? ','s')
if skip=='N'

shft=input('SHIFT THE CLIMATE SERIES FORWARD ONE YEAR? Y/N [N] ','s')
if isempty(shft), shft='N';   end

IA = E1(1,1):E1(1,2);    % Form some vectors of years.
IB = E1(2,1):E1(2,2);
IX = E1(3,1):E1(3,2);

if shft=='Y'
	IB=IB+1;  % Shift the climate series one year.
	IX(1)=[];  % Truncate the first year of analysis period.
else
end

[m1,n1] = size(A1);
[m2,n2] = size(B1);
m3=length(IX);   % number of years in analysis period

% Check for consistency of year coverage.

L1 = [length(IA)==m1 length(IB)==m2];
L2 = [(IX(1)>=IA(1) & IX(1)>= IB(1))   (IX(m3)<=E1(1,2) & IX(m3)...
	<=E1(2,2))];

if ~all(L1)
	error('LENGTHS OF A1,B1 DO NOT MATCH SPECIFIED YEARS IN E1');
end

if ~all(L2)
	error('SPECIFIED PD FOR CORRS NOT ENTIRELY COVERED BY A1, B1');
end


IAX = IA >= IX(1) & IA <= IX(m3);
IBX = IB >= IX(1) & IB <= IX(m3);

A1 = A1(IAX',:);  % Now A1 has only the specified years for analysis
B1 = B1(IBX',:);  % And so does B1.


else
end  % of loop on skip for quicky debug mode

%*****************   PREALLOCATE   ****************************


n4=length(I2);    % Number of valid tree sites for the analysis.
R=zeros(ndcols+1,ndcols+1);
R1=zeros (n4,4);
R1(:,1) = I2;   % Tree-ring vbles, relative to original array A1.
R2=zeros(n4,ndcols);  % 20 hard-coded to match ncols from near.m.

%*****************   COMPUTE CORRLN MTX   ***************


for i = 1:n4;  % loop for all tree-ring sites
	D2 = D1(I2(i),1:ndcols);  % Vector of station subscripts.
	ndnum=nD(I2(i));  % Number of neighbors for this site.

	DI = D2 <  zeros(1,length(D2));  % Index to invalid entries.

	D2(DI) = [];  % Remove  invalid members of the vector D2.
	
	aa=A1(:,I2(i));  % Pull proper tree-ring variable.

	R=corrcoef([aa B1(:,D1(I2(i),1:ndnum))]);  % corr mtx tree vs clim

	RSUB=	R(1,2:ndnum+1);  % corrs of this site with nearby stations.
	R1(i,2) = ndnum;  % Number of stns with corrlns for this site.
	R1(i,3) = median(RSUB);  % Median corrln for this site.
	R1(i,4) = max(RSUB);  % Maximum correl for this site.
	R2(i,1:ndnum) = RSUB;  % Store the correlns with nearest stns.
	
end

%   Make a surfer-usable xyz array   R4

if ~exist('tree')
	disp('LOAD tree.73, HOLDING THE X,Y COORDS')
	keyboard
else
end

R4=  [tree(I2,:)   R1(:,2:4)];
