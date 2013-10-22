function yx = lagyr1(X,Y,lg,yrs,xfill)
% Form a time series array of predictand and lagged values of predictors.
% Puts year in first column of output array, predictand in second column,
% unlagged predictors in next cols, lag -1 ,-2, etc, predictors, in next
% cols, lag +1, +2, ... predictors in next cols, and dummy values xfill
% in next cols. 

% D. Meko, 3-27-92

% lagyr1.m does not produce a long-term predictand/predictor matrix
% for use with reconstruction.  Usee lagyr2.m for that.

%******   INPUT ARGS
%
% X (mX X nX) unlagged array of predictors, possibly with a year column
% Y (mY x nY) array of predictands, possibly with a year column.
% lg(4 x 1)  shifting and lagging parameters
%    1 number of years to shift X relative to Y 
%	 2 number of negative lags to include on X
%    3 number of positive lags to include on X
%    4 number of desired variables per row in output array
% yrs (3 x 2)  beginning and ending years of 
%    row 1   Y
%    row 2   X
%    row 3   desired output array xy
% xfill (1 x 1)  dummy fill value to make output array rectangular

%*****   OUTPUT ARG
%
% yx (? x nprow) the array of Y and lagged X


%******   preallocate

myx=yrs(3,2)-yrs(3,1);  % # years in output array

%**********   Get rid of years columns of X, Y

k1=[]; k2=[];
k1=input('YEARS COLUMN IN X?  Y/N [N] ','s');
if isempty(k1), k1='N'; end
if k1=='Y', X(:,1)=[]; end
k2=input('YEARS COLUMN IN Y?  Y/N [N] ','s');
if isempty(k2), k2='N'; end
if k2=='Y', Y(:,1)=[];  end

[mX,nX]=size(X);  [mY,nY]=size(Y);


%****  Compute number of filler variables needed to make output 
%  time series array rectangular

nv1=  nX * (lg(2) + lg(3) + 1) + 2;
nfill=  lg(4)-rem (nv1,lg(4));
if nfill==lg(4), nfill=0; end
nv2=nv1+nfill;  %  total number of variables in output array

%*********  Make year pointers

yrY=(yrs(1,1):yrs(1,2))';
yrX=(yrs(2,1):yrs(2,2))';
yrYX=(yrs(3,1):yrs(3,2))';

LY= yrY >= yrs(3,1)  & yrY <= yrs(3,2);
L1=yrX >=  (yrs(3,1)+lg(1))  &  yrX <= (yrs(3,2)+lg(1));
LL=zeros(length(L1),(1+lg(2)+lg(3)));  % Initialize pointer to X
LL(:,1)=L1;  % first col points to unlagged rows of X

%  Build pointer to cols of X for negative lages
if lg(2) ~= 0 ;  %  if some negative lags are in model
	for i=1:lg(2);
		LSUB=L1;
		LSUB(1:i) = [];  % Lop off first i rows of pointer to X
		LL(:,i+1)=[LSUB;zeros(i,1)];
	end
end

% Append pointers to any positive lags on X
if lg(3) ~=0 ;  %  some positive lags on x are in model
	for i=1:lg(3);
		LSUB=L1;
		LSUB(length(L1)-i+1:length(L1))=[];
		LL(:,i+1+lg(2)) = [zeros(i,1); LSUB];
	end
end
size(LL)

xlagged=[yrYX  Y(LY,1)  X(L1,:)]; % First put year, predictand, and
%		unlagged X in first cols

% Append cols for lags
if (lg(2)+lg(3)) ~= 0
	for i=1:(lg(2)+lg(3));
		xlagged=[xlagged  X(LL(:,i+1),:)];
	end
end

% Append any additional cols needed to make matrix rectangular for any
% given year
if nfill~=0;  % If need dummy variables to get  rectangular array
	xlagged=[xlagged  xfill(ones(sum(L1),1),ones(nfill,1))];
else
end

% Rearrange output array to have desired number of cols
nn=nv2/lg(4);  % number of rows per case in target array
XL = zeros(lg(4), sum(L1)*nn);
size(XL)
XLT=(xlagged)';
XLS=XLT(:);
XL(:)=XLS;
yx=XL';
