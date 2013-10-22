% dwgt.m   inverse distance weighting of tree-ring network onto
%		a specified grid.  Assumes have input arrays with
%		coord locations for sites and grid points, and a data
% 		array of tree-rings at the sites.  Produces a new
% 		data array of tree-ring indices at the grid points.
%		Each row a year, each col a grid point.

%  	Modified 7-31-92 to allow restricting the number of sites to
% 		enter into the weighting.

%**************   PRELOADS  *****************************************

% F	? x 2, x-y coords of tree-ring sites
% G   ? x 2, x-y coords of grid points
% A   ? x mF, tree-ring data array at sites, each row a year

%************* OUTPUT ARRAY *****************************************

% Z   mA x mG, time series array of gridded tree rings

%******************* SIZE ARRAYS *********************************

[mF,nF]=size(F);
[mG,nG]=size(G);
[mA,nA]=size(A);

%***************  PREALLOCATE  **************************************

L = zeros(mG,1);
L1 = zeros(mG,1);
B = zeros(mG,mF);
mB=mG; nB=mF;

%******************** 	BEGIN  ****************************************

d= 213;   % hard-coded search distance
dcrit=30;  % minimim allowable radius for purpose of gridding.  Any
%		radius < dcrit  will be set to decrit before computing 
%		weights

d1=1/d;
d2=1/dcrit;

for i = 1:nA;    % loop for each tree-ring site

	K = F(i,:);  % Dupe a row of the coords of the tree-ring site
	L = K(ones(mG,1),:);  %  continue duping

	M = G - L;  % difference of x and y coords of site and each grid pt
	N = M .* M;  % square the delta x and delta y
	P = N(:,1) + N(:,2);  % sum the square of delta x and delta y
	Q = sqrt(P);  % the distance from the site to each of mG grid pts.
	U = 1 ./ Q;  % the inverse of the distance

	L2 = U >= d2(ones(mG,1),1); % if inverse of distance is greater
		% than 1/dcrit, set inverse distance to 1/dcrit.
	U(L2) = d2(ones(sum(L2),1),1);

	L1 = U >= d1(ones(mG,1),1);  % begin masking out distances > d.
	R = U .* L1;  % inverse distance, like U, but zero if < 1/d.


	B(:,i) = R;  % fill col of B with relevant inve dist for site i

end

%*******     LIMIT NUMBER OF STATIONS TO USE IN WEIGHTING?

nlim=input('LIMIT  NUMBER OF SITES WITHIN SEARCH DISTANCE? Y/N [N]',...
   's');
if isempty(nlim), nlim='N';  end
if nlim=='Y'
	maxnum=input('USE NEAREST THIS MANY SITES WITHIN SEARCH D:  '); 
	[mB,nB]=size(B);
	[Y,I]=sort(B');
	for i=1:mB
		B(i,I(1:nB-maxnum,i)) = zeros(1,nB-maxnum);
	end
else
end







S=(sum(B'))';  % col vectors of sums of invs distances for all sites for
	% a particular grid point.  If 93 grid points, 93 rows.

T= S(:,ones(nA,1));  % dupe the cv S into nA cols.

C = B ./ T;  % divide the row of distances for gpoint i by the sum
	% of the distances over all sites for that gridpoint
	% These are the inverse distance weights

Z = A * C';  % the distance weighted tree-ring data array
