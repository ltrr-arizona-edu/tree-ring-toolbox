% near1.m    by D. Meko, Dec 15, 1991

% Like near.m, but searches by quadrant.

% Given an input array of x,y coords of tree-ring sites, and a similar
% array of locations for climate stations, finds nearest n climate 
% stations in each of four quadrants to each tree-ring site, or finds
% nearest up to n stations within a specified search radius.  
% Also returns corresp distances.

% Assumes you have both types of data on the same x-y coordinate grid.
% Uses Pythagorean theorem for distance.  Adjusts from input units
% to other units (say, km) by an internal hard-coded factor--fscale.


%***************   PRELOADS  *********************************

% A (m1 x 2) [R] array of x,y coords for tree-ring sites
% B (m2 x 2) [R] likewise for climate data.


%*************8   SCREEN PROMPTED INPUT  ***************************

% opt (1 x 1) [L] 1=nearest n stations;  0=nearest "up to n" stations
%		not outside the search radius s.
% nqmax (1 x 1) [I] max number of stations to find for each quadrant
% ncols (1 x 1) [I]  number of desired columns in D1,D
% s (1 x 1) [R] search radius (assumed to be in real distance units)
% fscale (1 x 1) [R]  map scaling factor.  Say to convert from map
%		units to km.

%**************   RETURNED  *******************************************

% D1 (m1 x n3) [I] subscript array of nearest sites, corresp to cols of A 
% D (m1 x n3) [R] distances corresp to D1
% nD (m1 x 1) [I] number of valid stations for each tree site.
% I2 (m2 x 1) [I] pointer to valid gridpoints

%**************  OTHER KEY VARIABLES  ***********************************

% ncols - desired no. of cols in D1,D.  



%*******************   NOTES   ****************************************


% Usually save D1, D, nD and a cv I2 in a .mat file like TRP160.MAT,
% which would mean TRee centered to Precipi with 160-km search radius
 
%**********************  BEGIN  ***************************************

ABIN=input('A and B PRELOADED? Y/N [Y]','s');
if isempty(ABIN), ABIN='Y'; end;
if  ABIN=='Y'
else
	keyboard
end


s=input('SEARCH RADIUS:'  );
nqmax=input('MAXIMUM NUMBER OF STATIONS PER QUADRANT:  ');
opt=input('1=NEAREST n STATIONS;  0=WITHIN SEARCH DISTANCE ONLY:  ');
ncols=input('NUMBER OF DESIRED COLUMNS IN D1,D:  ');
clc,home;


a=-2.0;
ia=-2;

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
	disp(i)
	A1=ones(m2,1) * A(i,1:2);
	C1=A1-B(:,1:2);
	C2=C1 .* C1;
	C3=(sum(C2'))';
	C4=sqrt(C3);
	
	[Y,I]=sort(C4);
	C1R = C1(I,:);  %  Delta x, y , but ranked according to distance in C4

	Q1 = find(C1R(:,1) <0   & C1R(:,2) < 0 & Y(:,1)<=s);  % Quad 1 (upper right)
	Q2 = find(C1R(:,1) <0   & C1R(:,2) >= 0 & Y(:,1)<=s) ;
   Q3 = find(C1R(:,1) >=0  & C1R(:,2) >= 0 & Y(:,1)<=s);
   Q4 = find(C1R(:,1) >=0  & C1R(:,2) < 0 & Y(:,1)<=s);

	nq = [min([nqmax length(Q1)])   min([nqmax  length(Q2)])  ...
         min([nqmax length(Q3)])   min([nqmax  length(Q4)])];
	n = sum(nq);
   nD(i)=n;

	Y=Y * fscale;

	if nq(1)~=0, Q1=Q1(1:nq(1)); end;
	if nq(2)~=0, Q2=Q2(1:nq(2)); end
	if nq(3)~=0, Q3=Q3(1:nq(3)); end;
	if nq(4)~=0, Q4=Q4(1:nq(4)); end;

	if n>0
		D1(i,1:n) = [(I(Q1))'  (I(Q2))'  (I(Q3))'  (I(Q4))' ];
		D(i,1:n) = [(Y(Q1))'  (Y(Q2))'  (Y(Q3))'  (Y(Q4))' ];
	else
	end
	
	D1(i,n+1:ncols) = ia(ones(1,1),ones(ncols-n,1));
	D(i,n+1:ncols) = a(ones(1,1),ones(ncols-n,1));

	end
end


%if opt==0;  % Flag entries of D1 and D as invalid if outside search radius.
%	G=D>s;  % Logical=1 if stns outside search radius.
%	jj=sum(sum(G));  % Number of too-distant stations.
%	D1(G)=ia(ones(jj,1),:);
%	D(G)=a(ones(jj,1),:);
%	nD = nD - (sum(G'))';  % No. of stns within search distance.
%end

		
% This code added 11-18-92.  Compute cv I2, a pointer to rows of D, nD,D1
% specifying which gridpoints are valid for future gridding (by e.g., 
% dwgt1.m).  To be valid, must have found a station in the search radius.
 
I2=1:length(nD);  % initialize to assume all gridpoints valid
I2=I2(nD~=0);  % gets rid of rows with zero stns in search distance



