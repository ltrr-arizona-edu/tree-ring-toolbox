% dwgt1.m  inv-dist weighting using distance arrays previously calculated
%		by near1.m,near.m or near3.m.   Dwgt1.m differs greatly from dwgt.m.  
%		Dwgt.m  had to compute the distances from stations to gridpoints.
%		Dwgt1.m has those distances previously calculated.


%***************   PRELIM  *************************************

% Run either near.m or near1.m to get the arrays D, D1, I2, nD.
% These are usually stored in  .mat array, say trp160.mat -- indicating
% (tr)ee ring sites are "grid points", (p)recip stations are stations,
% and 160 km is the search radius.


%******************   INPUT   *************************************

% nD (m4,1) i  number of stations to be used in weighting each gridpoint.
%		Some elements might be 0, meaning cannot grid for this gridpoint.
%		m4  corresps to number of potential gridpoints.  In one applic, m4
%		is the number of tree-ring sites, around which representative ppt
%		series are to be weighted before single-site recon.   In another
%		applic, m4 is the number of 2 x 3 degree lat-long gridpoints to 
%		which station ppt is to be gridded.
%
% D1 (m3,n3) i  Which stations go into the weighting for each of
%		m3 gridpoints.  A row is for a potential gridpoint.  Cols are
%		pointers to cols of the data array A.  Rows are filled with
%		dummy values -2 for invalid stations.  So, for example, 
%		nD(3) = 5 means five stations go into gridding for a point,
%     and D1(3,:)=  343  761  12  34  45  -2 -2 -2 -2 -2 -2 -2
%		indicates which cols of A are the five stations.
%
% D (m2,n2) r  Distances from each gridpoint to each station.  Dummy
%		filled with -2.0.  D is same size as D1, and indicates corresp
%		distances for D1.  All rows of D,D1 may not be used in the gridding
%		depending on the index array I2.
%
% I2 (m1,1) i  pointer to cols of original array of potential gridpoints,
%		telling which are to be treated as valid in the analysis.  Points to
%		rows of D, D1, nD.
%
% A (m5,n5) r  The data to be gridded. m5 years by n5 stations.


%************   OUTPUT   *********************************************

% P (m5,m1) r  The gridded data.  m5 years, m1 gridpoints.


%**************   USE   **********************************************

% Convert the data array A to pctn (if ppt) or to deviations from mean
% (if temp) before the analysis.  These transformed variables are more
% coherent in space than the absolute values, and essentially remove 
% effects due to elevation.

% Best to run on not only the transformed data A, but on the 2-row array
% of station means and coeffs of variation.  If call the gridded
% means and cvs  PMN, PCV,  then store P, PMN, and PCV together in a .mat
% array.

% Can also apply to get nD-station mean (even weighting of stations) by
% specifying in response to prompt a dcrit larger than any distance in 
% D.  This action overrides the distances in D and treats all nD(i)
% stations for gridpoint i as if they are the same distance from the
% gridpoint.

% Right now, 4mb ram not enough to handle A 82x990 and W 990x244
% So, I have commented out last steps.  Do these separately after
% leaving the pgm. and clearing other garbage.


%*************   LOAD CHECK   **********************************

chk=input('A LOADED AND TRANFORMED? D,D1,I2,nD LOADED? Y/N [Y] ','s');
if isempty(chk), chk='Y';  end;
if chk ~= 'Y'
	keyboard
else
end


%******************   SIZE   *************************************

[m1,n1]=size(I2);
[m2,n2]=size(D);
[m3,n3]=size(D1);
[m4,n4]=size(nD);
[m5,n5]=size(A);
W=zeros(n5,m1);   % weights array



%***************   SCREEN INPUT   *******************************

dcrit=input('MIN DISTANCE FOR SENSITIVITY:  ');


%*********   HANDLE DISTANCES SHORTER THAN MIN FOR SENSITIVITY  ******

dc = dcrit(ones(m2,1),ones(n2,1));

DD=D;   % Dont want to change D in this pgm.
I3 = DD<dc;  % Pointer to els with dist shorter than dcrit.
sI3 = sum(sum(I3));  % Total num of elmts in D with distance shorter
		% than dcrit.  Note that dummy elements (-2.0) will also be
		% counted, but dont worry.
DD(I3)= dcrit(ones(sI3,1),:);  %  substitute min distance for elmts 
		% closer than dcrit to gridpoints.


%****************   LOOP FOR EACH VALID GRIDPOINT   ***************

for i = 1:m1;   % m1 is row size of I2
	disp(i)
	j=I2(i);  %  pointer to rows of D,D1, nD
	num = nD(j);  %  number of stations 
	W(D1(j,1:num),i) = (1 ./  DD(j,1:num))';
	w1=W(:,i);
	ws = sum(w1);
	W(:,i)=w1 / ws;
end


%P = A * W;
