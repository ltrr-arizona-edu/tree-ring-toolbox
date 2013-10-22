function [P,m,nv] = pctntsm(X,k1,YRS)
%
%  Transform a data array X into pctg-of-normal or 
%  departure-from-normal matrix, with control over period for
%	computation of means and over handling of NaNs
%
% Meko 8-5-96
%
%************************** IN ARGS *****************
%
% X (m1 x n1)r   data matrix, possibly with some NaNs
% k1 (1 x 2)i    options
%	k1(1) handling of means period
%		==1 means computed on valid data in all rows
%		==2 means computed on valid data as marked by rows
%	k2(2)  pct normal vs departures
%		==1 pct normal (dividsion)
%		==2 departure
% YRS (2 x 2)i  start, end years of 
%		row 1: matrix X
%		row 2: normalizing period for means
%
%************************** OUT ARGS ******************
%
% P (m1 x n1)r   pctg of normal or departure from normal matrix
%		corresponding to X.  See k1 for whether pctg or departures,
%		and whether all avail data or a subset of year used for
%		normals
%
% m (1 x n1)r  means used for normalizing X
% nv (1 x n1)i  sample sizes used for means
%
%
%**********************************************************


a=NaN;
[m1,n1]=size(X);
Y=a(ones(m1,1),ones(n1,1)); % Initialize output matrix

if m1~= (YRS(1,2)-YRS(1,1)+1);
	error('Row size of X inconsistent with YRS(1,:)')
end


if YRS(2,1)<YRS(1,1)  | YRS(2,2)>YRS(1,2)
	error('YRS(2,:) does not point to subset of YRS(1,:)')
end


%***************************
% Put data to be used for means in X1
yr1=(YRS(1,1):YRS(1,2))'; % year vector for X
LL = yr1>=YRS(2,1) & yr1<=YRS(2,2);  
if k1(1)==1; % use all rows for column means
	X1=X;
elseif k1(1)==2; % use a specified row subset for col means
	X1=X(LL,:);

else
end
%**************************

% Mark the NaN entries in X1, and replace these with zeros in 
% a copy-matrix X2
[m2,n2]=size(X1);
X2=X1;
L1=isnan(X1);
L2=~L1;
suma = sum(sum(L1));  % total number of NaNs in the rows period
if suma>0;
	X2(L1)=zeros(ones(suma,1),1);
end

%*****************************************************
% Compute column sums of X2
s2= sum(X2);



%***************************************************

% Compute number of valid data values (not NaN) in each col
% of the culled data matrix X2;  need this row vector to divide
% sums by to compute means
if suma==0;
	nv=m2(:,ones(n2,1)); % row vector
else
	nv =  sum(L2);
end

%*******************************************************

% Compute the means for the normalizing period
m = s2 ./ nv;

% Dupe the rv of means into full-sized matrix, same size as X 
M = m(ones(m1,1),:);

% compute pctg normal series or departures
if k1(2)==1; % division
	P=X ./ M;
elseif k1(2)==2;
	P=X - M;
else
	error('k(2) must be 1 or 2');
end
