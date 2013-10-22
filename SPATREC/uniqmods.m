function [L6,m]=uniqmods(Z,yrs)
% [L6,m]=uniqmods(Z,yrs)
%
% Pointer matrix to operative rows of tree-ring predictor
% matrix for ITI model
%
% D Meko  8-20-95
%
%
%*************** IN ARGS *************************
%
% Z (mZ x nZ)r tsm of tree indices
% yrs (3 x 2)i start and years in Z
%		row 1: entire Z
%		row 2: reconstruction period in Z
%		row 3: calibration period in Z
%
%
%*************  OUT ARGS ********************
%
% L6 (mL6 x nL6)L pointer matrix telling which cols of Z
%		are active for each of mL6 models
% m (mm x nm)i which model (row of L6) applies for each
%		of the reconstruction years in Z as pointed to by 
%		row 2 of yrs
%
%************ USER-WRITTEN FUNCTIONS CALLED
%
%
%************* NOTES *******************************
%
% For any year in the reconstruction period, find which trees are
% present and also cover all years in the calibration period
%
% Original application was to get a pointer index for a precip
% reconstruction for the San Pedro River Basin


[mZ,nZ]=size(Z);
yrall = (yrs(1,1):yrs(1,2))'; % cv of years of matrix Z


% Pointers to rows of calib and recon periods in Z
LZW = yrall >= yrs(3,1) & yrall <= yrs(3,2); % to calib period
LZX = yrall >= yrs(2,1) & yrall <= yrs(2,2); % to recon period

W = Z(LZW,:); % calibration subset
X = Z(LZX,:); % reconstruction subset
[mW,nW]=size(W);
[mX,nX]=size(X);

% Loop over each reconstruction year, finding the set of trees 
% (columns in X or Z) to be used as potential predictors
for k = 1:mX
	L1  = ~isnan(X(k,:)); % rv, trees present in year k
	% Exclude from this set any trees not covering the calib period
	W1 = W(:,L1); % subset of W -- only those key-year trees
	I1 = find(L1);
	L3 = isnan(W1); % matrix
	L4 = any(L3);  % rv, 0 if no values in the col NaN
	if any(L4);
		nz = sum(L4);
		L1(I1(L4)) = zeros(1,nz);
	end
	L5(k,:) = L1;
end

% Part one is now finished. We have L5 the logical pointer matrix
% telling which variables (cols of X) can be included as potential
% predictors for each reconstruction year.
%
% Each row of L5 indicates the potential predictors for a "model".
% Because the tree population does not change each year, many of
% these models are redundant. Now we will reduce the number of rows
% of L5 so that only unique models are included.  We will also build
% a col vector indicating which of these models applies in each
% reconstruction year.

% Store the indicator of models in the cv  m
% Initialize m to "1"
% 
L6=L5; % initialize matrix to hold reduced set of models
m = ones(mX,1);
mthis=1;
k3=m;
for i = 1:mX;
	if k3(i) ~=0;
		m(i) = mthis;
		LL1 = L5(i,:);
		LL1 = LL1(ones(mX,1),:);
		LL2 = LL1 == L5;
		LL3 =   (all(LL2'))'; % cv
		nn1 = sum(LL3); % number of years with this model
		k3(LL3) = zeros(nn1,1);
		m(LL3)= mthis(ones(nn1,1),:);
		L6(mthis,:)=LL1(1,:);
		mthis = mthis+1;
	end
end

L6 = L6(1:max(m),:);
