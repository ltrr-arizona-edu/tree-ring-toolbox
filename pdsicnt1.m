function [AC,BC,AI,BI,D,b]= pdsicnt1(A,B,yrs,a)

% Compares gridded PDSI with gridded tree rings to determine ability to
% identify droughts of specified severity.  Also tabulates "drought-area"
% series by counting number of gridpoints in drought each year.
% 
% You specify the threshold PDSI for drought.  For each gridpoint, this
% threshold value has an exceedance probability.  This function then 
% finds the tree-ring value for that gridpoint with the same exceedance
% probability, and uses that tree-ring value as threshold for tree-ring
% drought.

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
%		Calibration period is the overlap period you want to use for
%		comparing inferred with actual droughts.  Limited by possible
% 		nonexistence of data, and by poor quality, say in startup 4 yr of 
%		PDSI.
%
% a(1 x nA) threshold drought value for each of nA gridpoints.  
% 	 This would typically be a constant vector of some value.  Say 
% 	 a(1:nA)=-3.0.  Then only years with PDSI < -3.0 will classify as
%   drought years.


%**********  OUTPUT ARGS   **************

% AC (mA x 1)   PDSI count of number of points in drought each year
% BC (mB x 1)   tree count of number of points in drought each year
% D (5 x nA)   summary information on drought 
%   row 1:  number of years classified in drought- PDSI
%   row 2:  number of years classified in drought- trees
%   row 3:  number of 'hits' = droughts successfully diagnosed
%   row 4:  number of 'false droughts' = in trees but not PDSI
%   row 5:  number of  'missed' droughts = in PDSI, but not trees
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
D=zeros(5,nA);


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
%      INDEX CORRESP TO FIRST VALUE LARGER THAN CUTOFF THRESHOLDS
% 		 STORE SAME RANK VALUE FROM CORRESP COLS OF SORTED 
%	    TREE ARRAY BB.

for i = 1:nA
	I=find(AS(:,i)>a(i));
	b(i)= BS(I(1),i);
end


%**************   DUPE SOME VECTORS INTO MATRICES

aa=a(ones(mA,1),:);
bb=b(ones(mB,1),:);


%***************  COMPUTE INDICATOR ARRAYS AND TIME SERIES COUNTS

AI = A <= aa;
BI = B < bb;

AC = (sum((AI)'))';
BC = (sum((BI)'))';


%**********  HITS AND MISSES CHECK

AAI = AI(LA,:);
BBI=BI(LB,:)

D(1,:) = sum(AAI);  % number of drought yr in calib pd
D(2,:) = sum(BBI);
D(3,:) = sum(AAI & BBI); % number of hits
D(4,:) = sum (BBI  &  (~AAI));  % false alarms
D(5,:) = sum (AAI &  (~BBI));  % missed droughts
D(6,:) = D(3,:) - D(5,:);   % hits minus misses
