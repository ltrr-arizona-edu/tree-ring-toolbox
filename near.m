% near.m    by D. Meko, Sept 24, 1991

% Given an input array of x,y coords of tree-ring sites, and a similar
% array of locations for climate stations, finds nearest n climate 
% stations to each tree-ring site, or finds nearest up to n stations
% within a specified search radius.  Also returns corresp distances.

% Assumes you have both types of data on the same x-y coordinate grid.
% Uses Pythagorean theorem for distance.  Adjusts from input units
% to other units (say, km) by an internal hard-coded factor--fscale.


%***************   PRELOADS  *********************************

% A (m1 x 2) [R] array of x,y coords for tree-ring sites
% B (m2 x 2) [R] likewise for climate data.


%*************8   SCREEN PROMPTED INPUT  ***************************

% opt (1 x 1) [L] 1=nearest n stations;  0=nearest "up to n" stations
%		not outside the search radius s.
% n (1 x 1) [I] number of stations to find for each station
% ncols (1 x 1) [I]  number of desired columns in D1,D
% s (1 x 1) [R] search radius (assumed to be in real distance units)
% fscale (1 x 1) [R]  map scaling factor.  Say to convert from map
%		units to km.

%**************   RETURNED  *******************************************

% D1 (m1 x n3) [I] subscript array of nearest sites, corresp to cols of A 
% D (m1 x n3) [R] distances corresp to D1
% nD (m1 x 1) [I] number of valid stations for each tree site.
% I2 (m2 x 1) [I] pointer to valid gridpoints (those with nonzero nD)

%**************  OTHER KEY VARIABLES  ***********************************

% ncols - desired no. of cols in D1,D.  



%*******************   NOTES   ****************************************

 
%**********************  BEGIN  ***************************************

ABIN=input('HAVE YOU LOADED A and B? Y/N [Y]','s');
if ~isempty(ABIN), ABIN='Y'; end;
if ABIN=='N'
	keyboard
else
end


s=input('SEARCH RADIUS:'  );
n=input('MAXIMUM NUMBER OF NEARBY SITES TO CONSIDER:  ');
opt=input('1=NEAREST n STATIONS;  0=WITHIN SEARCH DISTANCE ONLY:  ');
ncols=input('NUMBER OF DESIRED COLUMNS IN D1,D:  ');
clc,home;


a=-9.0;
ia=-9;

%*********   CHANGE FROM DEFAULT SCALING OF UNITS  ****************

fscale=1.000;
disp(['CURRENT FSCALE = ',num2str(fscale)]);
fch=input('CHANGE FSCALE? Y/N [N]','s');
if isempty(fch), fch='N';  end;
if fch=='Y'
	fscale=input('NEW VALUE FOR FSCALE:  ');
else
end


[m1,n1]=size(A);
[m2,n2]=size(B);

D1=zeros(m1,ncols);
D=zeros(m1,ncols);

for i=1:m1

	A1=ones(m2,1) * A(i,1:2);
	C1=A1-B(:,1:2);
	C2=C1 .* C1;
	C3=(sum(C2'))';
	C4=sqrt(C3);
	
	[Y,I]=sort(C4);
	Y=Y * fscale;

	D1(i,1:n) = I(1:n)';
	D(i,1:n) = Y(1:n)';


	end
end

ndiff=ncols-n;
nD=ones(m1,1) * n;
D1(:,n+1:ncols) = ia(ones(m1,1),ones(ndiff,1));
D(:,n+1:ncols) = a(ones(m1,1),ones(ndiff,1));


if opt==0;  % Flag entries of D1 and D as invalid if outside search radius.
	G=D>s;  % Logical=1 if stns outside search radius.
	jj=sum(sum(G));  % Number of too-distant stations.
	D1(G)=ia(ones(jj,1),:);
	D(G)=a(ones(jj,1),:);
	nD = nD - (sum(G'))';  % No. of stns within search distance.
end


% This code added 11-18-92.  Compute cv I2, a pointer to rows of D, nD,D1
% specifying which gridpoints are valid for future gridding (by e.g., 
% dwgt1.m).  To be valid, must have found a station in the search radius.
 
I2=(1:length(nD))';  % initialize to assume all gridpoints valid
I2=I2(nD~=0);  % gets rid of rows with zero stns in search distance


