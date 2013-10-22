function [P,Z]=pairchk1(X);
% Proportion of daily ppt values missing or zero for various Julian days
% and stations -- poorly vectorized version;  see pairchk2.m for fully
% vectorized
%
% X (mX x nX) matrix of daily ppt at several stations.  Col 1 assumed to be
%	year,  col 2 the Julian day, remaining cols the ppt for several stns
%
%
% P (366 x nX-2) number of non-missing values for each Julian day at each stn
% Z (366 x nX-2) proportion of non-missing values that are zero
%
%
%******* D. Meko   11-1-93
%
%
%*****  CONTEXT  ************
%
% One of many functions under top-level function dayfill.m, used to estimate
% missing daily ppt in a matrix of ppt data
%
%
%*********  GLOBALS -- none
%
%

[mX,nX]=size(X);

P=zeros(366,nX-2);
Z=zeros(366,nX-2);


% Loop for each Station
for i1 = 3:nX;
	i2=i1-2;
	disp(['Now on station ',int2str(i2)]);
	xday=X(:,2);
	x=X(:,i1); % pull data for one station
	% Loop for each Julian day
	for j=1:366;
		disp(['Now on Julian day ',int2str(j)]);
		i3=xday==j;  % logical vector to rows for this Julian day
		x1=x(i3);   % Data for this Julian day
		i4=x1>=0;  % logical pointer to rows with non-missing data
		P(j,i2)=sum(i4);
		x2=x1(i4);  % subset vector of non-missing data
		i5=x2==0;  % logical pointer to zero values in subset vector
		Z(j,i2)=sum(i5)/ sum(i4);
	end;  % of loop for julian day
end;  % of loop for stations


