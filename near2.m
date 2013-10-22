% near2.m  Find nearest climate stations to gridpoints.

% Given arrays of x,y coords for gridpoints and stations, finds out which
% gridpoints are nearest each station.  Forms arrays pointing to which 
% stations, and the distances, for each gridpoint.  

% Assumes gridpoint and station coords are in same x,y system, and gets
% distance by Pythagorean theorem.  So be sure have an equal-area
% projection.  I use DISPLAA on convex to convert long-lats to map units, 
% then map2km.f to convert map units to km before running this program.

% First used to form distance arrays for later gridding of hcn precip
% data.


%******************   PRELOADS   **************************

% A (m1 x 2)  x,y coords of gridpoints
% B (m2 x 2)  x,y coords of stations


%***************   SCREEN PROMPTS   ************************

% fscale (1 x 1) r  map-scaling factor, default value 1.0.  Just in case
	%  you don't have x,y coords in desired units.
% ncols (1 x 1) i  number of desired cols in output arrays D, D1.  Make at
	% least as large as largest # valid stations for a
	% gridpoint.
% opt (1 x 1) L  option array not yet needed
% s (1 x 1) r  max allowable distance, gridpoint to station


%***************   OUTPUT   **********************************

% D1 (m1 x ncols) i  Pointer to B telling which stations go with this
	% gridpoint.
% D (m1 x n3) r  Corresponding distances for D1.
% nD (m1 x 1) i  Valid number of stations for each gridpoint.



%***************   LOAD CHECK   *****************************

ABIN = input('A and B pre-loaded? Y/N [Y] ','s');
if ~isempty(ABIN), ABIN='Y';  end;
if ABIN == 'N'
	keyboard
else
end


%**************   SCREEN INPUT  *********

s=input('MAX ALLOWABLE DISTANCE, SITE TO GRIDPOINT: ');
ncols = input('# COLS IN D1, D:  ');  % Make this at least as
	% large as the largest number of stations that will
	% fall within s of any gridpoint.  May need a little 
	% trial and error.


%*************   OPTIONAL SCALING OF DISTANCE   ******************

fscale = 1.0;   % hard-coded
disp(['CURRENT SCALE = ',num2str(fscale)]);
fch = input('CHANGE SCALE? Y/N [N]','s');
if isempty(fch), fch='N'; end;
if fch=='Y'
	fscale=input('NEW VALUE FOR FSCALE: ');
else
end


%************   SIZE AND PREALLOCATE   ****************************

[m1,n1] = size (A);  % Size of gridpoint coord array.  m1=number grdpts
[m2,n2] = size(B);  % Size of station coord array.
a = -2.0;  % hard-coded dummy value to fill invalid cells in D.
ia = -2;  % hard-coded dummy value to fill invalid cells in D1.
D1 = ia(ones(m1,1),ones(ncols,1));
D = a(ones(m1,1),ones(ncols,1));
nD = zeros(m1,1);


% Treat each station in a loop.  Find distance to each gridpoint, then
% sort to find nearest gridpoint.  Then increment the counter nD for
% that gridpoint, and add the approp entries to D and D1.  

for i = 1:m2;  % For each station
	disp(i)
	B1= ones(m1,1) * B(i,1:2);  % Dupe coords for station into mult rows
	C1 = B1-A(:,1:2);  % Delta x, delta y to each gridpoint.
	C2 = C1 .* C1;  % Square of delta x, delta y
	C3 = (sum(C2'))';  % sos of delta x, delta y
	C4 = sqrt(C3);  % Corresp distances
	[Y,I] = sort(C4);  % Rank distances to gridpoints; I(1) holds index
		% of nearest gridpoint.
	I1= I(1);
	Y = Y * fscale;
	nD(I1)=nD(I1)+1;
	D1(I1,nD(I1))=i;
	D(I1,nD(I1)) = Y(1);
end


% Treat each gridpoint in a loop.  Sort D,D1 nearest station to 
% farthest, for all valid stations.

%for i = 1:m1;   
%	[y,isub] = sort(D(i,1:nD(i)));
%	D(i,1:nD(i)) = y;
%	D1(i,1:nD(i)) = D1(i,isub);
%end


% Set cells of D,D1 corresp to stations outside the max allowable 
% distance to the dummy values a, ia.

%G = D>s;   % Outside max allowable distance range
%jj = sum(sum(G));  % Num of elements of D outside ranges.
%D1(G) = ia(ones(jj,1),:);
%D(G) =   a(ones(jj,1),:);
%nD  = nD - (sum(G'))';

