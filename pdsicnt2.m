function [AC,BC,AI,BI,D,b]= pdsicnt1(A,B,yrs,a)

% Tabulates number of gridpoints with tri below a threshold.  Threshold
% is based on analysis of gridpoint pdsi.  For given pdsi value (e.g.,
%-2), number of drought years in  a "calibration period" is counted for
% each gridpoint.  Median number of drought years is computed over
% gridpoints.  This scalar number is then taken as the number of 
% tree-ring drought yrs at each gridpoint.  Threshold level for each grid
% point is set to give this number of tree-ring drought years.

% Also tabulates "drought-area" series by counting number of gridpoints 
% in drought each year.
% 
%*******   INPUT ARGS
%
% A (mA x nA)  PDSI, mA years and nA gridpoints
% B (mB x nB)  tree-ring index, mB years and nB gridpoints.  
%	  		   nA must = nB
% yrs(1 x 4)  1. first year of A
%             2. first year of B
%             3. first year of "calibration" period
%             4. last year of "calibration" period
%
%		Calibration period is the overlap period on which the pdsi data 
%     will be analyzed for frequency of drought. 
%
% a(1 x nA) threshold drought value for each of nA gridpoints.  
% 	 This would typically be a constant vector of some value.  Say 
% 	 a(1:nA)=-3.0.  Then only years with PDSI < -3.0 will classify as
%   drought years.


%**********  OUTPUT ARGS   **************

% AC (mA x 1)   PDSI count of number of points in drought each year
% BC (mB x 1)   tree count of number of points in drought each year
% D (2 x nA)   summary information on drought for each gridpoint
%   row 1:  number of years classified in drought- PDSI
%   row 2:  number of years classified in drought- trees
% AI (mA x nA) 1/0 indicator array for drought in PDSI array A
% BI (mB x nB) 1/0 indicator array for drought in tree array B
% b  (1 x nB)  threshold for tree drought for each gridpoint.  Growth
%    less than this is judged a tree drought

%************  SIZE SOME ARRAYS

[mA,nA]=size(A);
[mB,nB]=size(B);


%************  PREALLOCATE

b=zeros(1,nB);
AC=zeros(mA,1);   BC=zeros(mB,1);
AI=zeros(mA,nA);  BI=zeros(mB,nB);
D=zeros(2,nA);


%*********  COMPUTE YEAR VECTORS FOR PLOTS

yr1=(yrs(1):mA+yrs(1)-1)';     % for A
yr2=(yrs(2):mB+yrs(2)-1)' ;    % for B
yr3=  (yrs(3):yrs(4))';       % calibration period


%********  COMPUTE 0/1 VECTORS INTO A, B, POINTING TO CALIB PD YEARS

LA = yr1 >= yrs(3)  & yr1 <= yrs(4);
LB = yr2 >= yrs(3)   & yr2 <= yrs(4);



%*********  PULL CALIB SUBSETS FROM A, B, AND SORT SUBSET OF A

AA = A(LA,:);
BB = B(LB,:);
AS = sort(AA);
BS = sort(BB);


%**********  LOOP THRU COLS OF SORTED PDSI SUBSET ARRAY TO FIND
%      NUMBER OF YEARS IN PDSI DROUGHT AT EACH GRIDPOINT.

for i = 1:nA
	I=find(AS(:,i)>a(i));  %I(1)-1 yrs will have had a pdsi value
  	% less than or equal to a(i)  at this gridpoint
	ndrt(i)=I(1);  % number of drought years plus 1 
end
nmed=median(ndrt) % median of (number of years +1) 


%*************  TABULATE TREE-RING THRESHOLD FOR EACH GRIDPOINT.
%		THRESHOLD CORRESPONDING TO nmed YEARS IN DROUGHT IN CALIB PD.

for i=1:nB
	b(i)=BS(nmed,i);  % tri lower than b(i) means drought
end
	

%**************   DUPE DROUGHT-THRESHOLD VECTORS INTO MATRICES

aa=a(ones(mA,1),:);
bb=b(ones(mB,1),:);

%******  COMPUTE INDICATOR ARRAYS AND TIME SERIES COUNTS
%        FOR FULL LENGTH OF TIME SERIES (NOT JUST CALIB PERIOD)    
AI = A <= aa;  
BI = B < bb;

AC = (sum((AI)'))';
BC = (sum((BI)'))';


%**********  FILL SUMMARY ARRAYS

AAI = AI(LA,:);   % indicator of pdsi drought years, calib period
BBI=BI(LB,:);      % "            tri

D(1,:) = sum(AAI);  % number of drought yr in calib pd
D(2,:) = sum(BBI);
