function W=merg1(U,V,tailyrs)
% Merge an earthinf ppt file with tail end of another
% Written so that TRL last few years could be spliced on a 
% long record derived from earthinfo CD
%
%
%********************* IN ARGS
%
% U (mU x nU)r   earthinfo file (year and 12 monthly values in inches) 
%			Missing values as -9.99
% V (mV x nV)r   trl file (year and 12 or 13 values in hundredths
%		of inches;  13 values if annual total included
%		Missing values as -999
% tailyrs (1 x 2)i  start, end year of the part of TRL data to merge 
%		as end of W
%
%
%
%******************* OUT ARGS
%
%  W (mW x nW)r  the merged file
%
%
%***********************888  NOTES ************
%
% Make sure the segment of the TRL file represented by V has no 
% blank values.  In other words, manually substitute  -999s




[mU,nU]=size(U);
[mV,nV]=size(V);

% Lopp off annual total, if necessary, from TRL file
if nV > 13
	V=V(:,1:13);
end


% Pointers to years
LU1= U(:,1)<tailyrs(1);  % pointer to rows of U to use
LV1= V(:,1)>= tailyrs(1); % pointer to rows of V to use

% Pull out the segment of TRL data to use, and convert to inches
Uadd = V(LV1,:);
Uadd (:,2:13) = Uadd(:,2:13)/100;

% Build merged file
W = [U(LU1,:);  Uadd];


% Check that "year" column correct
yr = W(:,1);
yrdiff = diff(yr); % first difference
L2= yrdiff == 1;  % 1 year increment in years column
if ~all(L2)
	error('The year column is not sequential')
end
