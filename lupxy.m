% lupxy.m  :  build lookup table for x,y map coords from long, lat
%
% D. Meko 5-11-94
%
% Automates procedure I used to do stepwise.  Problem is how to plot
% sites on a .dwg map that has its own x,y coordinate system.  Sites
% have long, lat.  Solution is as follows:
%
% 1. digitize the 1:1-million map, usually after cutting up into parts
% 2. using dxfout, output separate files with lines from each degree to
%    next into separate dxf files for each part
% 3. using dxfrd2.f, convert the dxf files to .dat files:  2cols, 
%    corresp to x coord, y coord
% 4. Those .dat files are used by lupxy.m
%
%
% ************  Steps before running lupxy  ***************************
%
% 1. rename data files as Z1,Z2,Z3,....  Start at lower left of
%    study area (SW), work E across southernmost files, then work
%    W to E across next sub-maps to north, etc.
% 2. Make nZ x 4 matrix of beg and end long, lat for each of nZ
%    data files.  cols are 1=beg lon, 2=end long, 3=beg lat, 4=end lat
%
%
% *********************** script lupxy.m does this *****************
%
% 1. From info on bounding longs, lats, computes numb of rows cols each
%    submatrix (reshaped original data files)
% 2. Splits the nZ matrices into nZ single-col vectors of x and y;
%	  reshapes these into matrix blocks that will fit into tblx, tbly.
%    components
% 3. Initializes the tblx and tbly output lookup tables
% 4  Within a loop for each part:
%		- load a matrix representing a block of long-lat
%		- reshape the X or Y matrix,
%		- Put the X or Y matrix into its correct slot in tblx or tbly
%
%
%************************  NOTES  ********************
%
% Hard-coded initial section must be changed for each application
%



%******************* HARD CODE *******************************

nZ=12 ; % number of separate parts of map (x-y data files dxf'd)

%  begin, end longitude, latitude of each part
B=[120 116 28 32; 115 112 28 32;  111 109 28 32;  108 105 28 32;...
   120 116 33 35; 115 112 33 35;  111 109 33 35;  108 105 33 35;...
   120 116 36 38; 115 112 36 38;  111 109 36 38;  108 105 36 38]
B(:,1:2)=-1*B(:,1:2);  % make longitude negative so axis increases
%   to right


%Ranges of row and col indices into tblx for placing reshaped parts
PR= [2 6; 7 10; 10 13;  14 17];
PR=[PR;PR;PR];  % row index into tblx of each sub-matrix
PC = [2 6; 7 9; 10 12];
PC=PC(ones(4,1),:);
PC=reshape(PC,12,1); % column index
  

% Boundaries of the study area
lowlat=28;
hilat=38;
lowlon=-120;
hilon=-105;

%***********************  END HARD CODE

% Compute number of rows, cols in each part in what will be the reshaped
% data matrices
nrX=B(:,2)-B(:,1)+1;
ncX=B(:,4)-B(:,3)+1;



% initialize tblx and tbly
mT=lowlon-hilon+2; % number of rows for tblx, tbly
nT=lowlat-hilat+2; % number of cols for tblx, tbly

tblx=zeros(mT,nT);
tbly=zeros(mT,nT);
tblx(:,1)=[lowlon lowlon:hilon];
tbly(1,2:nT) = [lowlat hilat];



% Load the data matrices and do the work
% Read and process files with names Z1.dat, Z2.dat, ...
for k=1:nZ
	Zk = ['Z' int2str(k)];
	filename = [Zk '.dat'];
	if ~exist(filename), break, end
	eval(['load ' filename]);
	Z=eval(Zk);
	X=Z(:,2);
	Y=Z(:,1);
	nc=ncX(k);
	nr=nrX(k);
	Y=reshape(Y,nr,nc);
	X=reshape(X,nr,nc);
	tblx(PR(k,1):PR(k,2),PC(k,1):PC(k,2))=X;
   tbly(PR(k,1):PR(k,2),PC(k,1):PC(k,2))=Y;
end

