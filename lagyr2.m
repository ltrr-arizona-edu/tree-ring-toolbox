function [yx,yxf] = lagyr2(X,Y,lg,yrs,xfill)
% Make predictand/predictor arrays for use by BMDP
% Makes time series arrays of predictand and lagged values of predictors.
% Puts year in first column, predictand (or dummy value) in second column,
% unlagged predictors in next cols, lag -1 ,-2, etc, predictors, in next
% cols, lag +1, +2, etc predictors in next cols, and dummy values xfill
% in next cols. 

% D. Meko, 4-5-92

%*****   USAGE
%
% Typically, call function, save yx as ascii, and
% save yxf in a mat file for later use in generating reconstructions.
% Then run BMDP to get regression coefficients, import the coefs into
% matlab, and apply to the array yxf to get reconstruction.

%******   INPUT ARGS
%
% X (mX X nX) unlagged array of predictors, possibly with a year column
% Y (mY x 1 or 2) matrix predictand, possibly with a year column.
% lg(4 x 1)  shifting and lagging parameters
%    1 number of years to shift X relative to Y  (usually set = 0)
%        If nonzero, lg(1) should be negative, meaning tree for year
%        t+1 offset against climate for year t
%	  2 number of negative lags to include on X
%    3 number of positive lags to include on X
%    4 number of desired variables per row in output array.  This
%		   allows flexibility when large number of variables would
%        exceed 255 character record length if strung out in a single
%	      row.
% yrs (3 x 2)  beginning and ending years of 
%    row 1   Y
%    row 2   X
%    row 3   desired output array xy
% xfill (1 x 1)  dummy fill value to make output array rectangular

%*****   OUTPUT ARG
%
% yx (? x nprow) the array of Y and lagged X rearranged in BMDP form.
%		Covers only period given by yrs(3,:).
% yxf (? x ?) full-period version of yx, not in BMDP form.
%   Preditand/predictor matrix for period of available reconstruction
%   given matrix X and the lag and shift parameters.  yxf is multiplied
%   by regression coefs to get long-term reconstruction.


%**********   Get rid of years columns of X, Y

k1=[]; k2=[];
k1=input('YEARS COLUMN IN X?  Y/N [N] ','s');
if isempty(k1), k1='N'; end
if k1=='Y', X(:,1)=[]; end
k2=input('YEARS COLUMN IN Y?  Y/N [N] ','s');
if isempty(k2), k2='N'; end
if k2=='Y', Y(:,1)=[];  end

[mX,nX]=size(X);  [mY,nY]=size(Y);
if nY~=1, error('lagyr1.m handles only one predictor.'); end;

% Compute valid year coverage of reconstruction, from coverage by X and
% shift and lag parameters

yrfon=max([yrs(2,1)-lg(1)  yrs(2,1)+lg(1)+lg(2)]);  % X coverage and lag/shift
yrfoff=min([yrs(2,2)+lg(1)  yrs(2,2)+lg(1)-lg(3)]); % allow this recon yr cover

%yrfon= yrs(2,1)+lg(1)+lg(2);
%yrfoff= yrs(2,2)+lg(1)-lg(3); 

%****  Compute number of filler variables needed to make output 
%  time series array rectangular

nv1=  nX * (lg(2) + lg(3) + 1) + 2;
nfill=  lg(4)-rem (nv1,lg(4));
if nfill==lg(4), nfill=0; end
nv2=nv1+nfill;  %  total number of variables in output array

%*********  Make year pointers

yrY=(yrs(1,1):yrs(1,2))';
yrX=(yrs(2,1):yrs(2,2))';
yrYX=(yrs(3,1):min([yrfoff yrs(3,2)]))';
yrall=(yrfon:yrfoff)';

LY= yrY >= yrs(3,1)  & yrY <= max(yrYX);

L1=yrX >=  yrfon-lg(1)  &  yrX <= yrfoff-lg(1);
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

yrpoint= yrall >= yrs(3,1)  & yrall <= max(yrYX);
numyr=sum(yrpoint);
YY=xfill(ones(length(yrall),1),:);
YY(yrpoint)=Y(LY,1);


yxf=[yrall  YY  X(L1,:)]; % First put year, predictand, and
%		unlagged X in first cols

% Append cols for lags
if (lg(2)+lg(3)) ~= 0
	for i=1:(lg(2)+lg(3));
		yxf=[yxf  X(LL(:,i+1),:)];
	end
end
% Append any additional cols needed to make matrix rectangular for any
% given year
if nfill~=0;  % If need dummy variables to get  rectangular array
	yxf=[yxf  xfill(ones(sum(L1),1),ones(nfill,1))];
else
end


% Form a subset array in BMDP-compatible format.  Subset typically covers 
% the years for calibration/verification, as specified by yrs(3,:), 
% adjusted for possible loss of end years if no tree-ring coverage .

% First point to subset of years, then 
% rearrange output array to have desired number of cols per record.

LL3=yxf(:,1) >= yrs(3,1) & yxf(:,1) <= yrs(3,2); 
nn=nv2/lg(4);  % number of rows per case in target array
short=yxf(LL3,:);

XL = zeros(lg(4), sum(LL3)*nn);
XLT=(short)';
XLS=XLT(:);
XL(:)=XLS;
yx=XL';
