% Spatial reconstruction using no lags in single-site  models.
% Template was spatrec1.m.  Reconstruct values of a climate variable at
% regions or gridpoints from tree-ring indices of individual trees
%
% 
%******************* IN ARGS *****************************
%
% T (n x nt)r  tree indices, n years at nt trees; NaN filled
% C (mc x nc)r climate for mc years at nc stations; NaN filled
% ICD (m1 x n1)i  index matrix telling which stations are to averaged
%		to form the site-centered climate series.  m1==nt
%		For example, if climate stations 3,6 and 9 apply to tree #1,
%		row 1 of ICD would be:  3  6  9  NaN  NaN ...
%		The col size of ICD must be large enough to accomodate the 
%		largest number of climate stations assigned to any one tree.
%		For example, the single-site reconstruction would be of PPT
%		averaged over stations 3,6 and 9 from the tree index at
%		tree #1.
% ICY (m3 x n3)i for each region (row), which stations to average over
%		in forming regional climate series. NaN filled.  m3 is the
%		desired number of regions.  The number of valid (not NaN) 
%		entries in a row is the number of stations to be used for
%		the regional average
% IDX (m4 x n4)i for each predictor location (row), which single-tree
%		reconstructed variables to average together. "Which" refers to
%		column of T (or equivalently of D) 
%		group
% YRS (3 x 2)i beginning and ending years of:
%		row 1: T, the tree-index matrix
%		row 2: C, the station-climate index
%		row 3: the calibration period for the model;
%		row 4: reconstruction period: should not overlap calib period
% mineig (1 x 1)i minimum number of predictor amplitudes to use
%		If model has this or fewer predictor variables available,
%		all predictor eigs are included as potential predictors,
%		not just the eigenvalue-1 ones
% O (1 x nO):  options vector
%		O(1): AR modeling of tree-indices
%			0 = none;  1 = ARMA(1,1);   2=lowest AIC AR model
%
%
%************************  OUT ARGS ***************************
%
% Z (mZ x nZ)r  single-site reconstructions
% Y (mY x nY)r  regional reconstructions
% SZ -- stats on single-site reconstructions
% SY -- stats on regional reconstructions%
%
%
%
%************************  NOTES ********************************
%
% Any desired ARMA "whitening" of tree-ring series is assumed to
% already have been done before calling this function.  A typical
% approach might be to use the "residual" tree indices as output
% by ARSTAN or an equivalent program
%
% Station-climate series is assumed not to have NaNs in data
% used by any of the tree-centered or regional climate series
%
% standard error of prediction for regional recons
%
%
%
%******************* RATIONALE  *******************************
%
% Two stages of regression in attempt to handle the problem of
% variable autocorrelation of tree-ring series in a spatial-recon
% context. 
%
% Tree-based reconstruction rather than site-based because 
%		1. Some tree grouping other than that used in forming
%			site chronologies might be more appropriate for 
%			capturing the climate signal in the region.  For
%			example, several runoff-sensitive trees might be
%			scattered over several sites in a river basin.  Why
%			dilute their hydrologic signal by averaging together
%			with other trees at the sites?
%		2. Allows recon to extend back to the start year of the
%			oldest tree;  reconstructions take advantage of
%			whatever trees available in any year; validation
%			allows estimation of prediction accuracy, which will
%			vary year by year


%************************ OUTLINE OF PROGRAM STEPS **************
%
% 1 Check input arguments, size and allocate variables
% 2 Make site-average and regional average climate series
% 3 Single site reconstructions
% 4 Spatial reconstructions
%     * calib step
%			- find col-index subsets for all recon years
% 5 Crossvalidation
%		- single-tree regression
%		- spatial regression

%************************* DERIVED VARIABLES
% 
% ns1 -- cv, number of clim stations averaged around each tree
% ns2 -- cv, number of clim stations averaged in each region



%****************** STEP-1: ARGUMENT CHECK, SIZING, ALLOCATION


kdope=10;
if kdope~=1;

[n,nt]=size(T);
[mc,nc]=size(C);
[m2,n2]=size(YRS);
[m1,n1]=size(ICD);
[m3,n3]=size(ICY);
[m4,n4]=size(IDX);

if m2~=4| n2~=2, error('YRS must be 4 x 2'), end
if nt ~= m1, error('col size of T must equal row size of ICD'), end


if YRS(3,2) > YRS(1,2),
	error('End year of calib pd later than end year of tree rings')
end
if YRS(1,1)>YRS(3,1);
	error('Start year of T too late for specified calib pd')
end
if YRS(1,2)<YRS(3,2),
	error('End year of T too early for specified calib pd')
end


if YRS(2,1)>YRS(3,1)
	error('Start year of C too late for specified calib pd')
end
if YRS(2,2)<YRS(3,2),
	error('End year of C too early for specified calib pd')
end

if  YRS(4,1)<YRS(1,1),
	error('Recon period too early first year')
end
if YRS(4,2)>= YRS(3,1),
	error('Recon period overlaps calibration period')
end

% Some key year-index pointers
yr1 = (YRS(1,1):YRS(1,2))'; % years vector for T
yr2 = (YRS(2,1):YRS(2,2))'; % years vector for C
yr3 = (YRS(3,1):YRS(3,2))'; % years vector for calib period
yr4 = (YRS(4,1):YRS(4,2))'; % years vector for recon period

L1T = yr1>=YRS(3,1) & yr1<=YRS(3,2); % to cal pd years in T
L1C = yr2>=YRS(3,1) & yr2<=YRS(3,2); % to cal pd years in C

L2T = yr1>=YRS(4,1) & yr1<=YRS(4,2); % to recon pd years in T


%****************** STEP-2: MAKE SITE-AVE AND REGIONAL-AVE CLIM SERIES


% AVERAGE CLIMATE SERIES AROUND TREE SITES
% Number of stations to average around each tree
L1=~isnan(ICD);
if n1==1; % special case, ICD is not a matrix but a cv
	ns1=L1;
else
	ns1=   (sum((L1)'))'; % col vector of number of clim station each tree
end
if any(ns1<1), error('Some row of ICD is all NaN'), end

% Form logical matrix corresponding to ICD that can be used in 
% matrix multiplication 
L2 = ind2log(ICD,nc);
D = C * L2';  % tsm of sum of ppt over stations near each site

% Check for NaNs in selected series
if any(any(isnan(D(L1C,:)))),
	error('NaN found in calib part of single-tree climate matrix D')
end

temp1=ns1'; % convert to row vector -- number of stations each site
denom=temp1(ones(mc,1),:); % expand to matrix
D = D  ./ denom;  % climate variable averaged over specified
		% stations near each tree;  D is the tree-centered climate
		% series to be used as the predictand variable in single-
		% tree reconstructions



% AVERAGE THE CLIMATE SERIES INTO "REGIONS" OR GRIDPOINTS TO BE
% RECONSTRUCTED
% Number of stations to average in each region
L1=~isnan(ICY);
if n3 == 1; % special case: ICY is not a matrix but a cv
	ns2=L1;
else
	ns2=   (sum((L1)'))'; % col vector of number of clim station each tree
end
if any(ns2<1), error('Some row of ICY is all NaN'), end

% Form logical matrix corresponding to ICY that can be used in 
% matrix multiplication 
L2 = ind2log(ICY,nc);
Y = C * L2';  % tsm of sum of ppt over stations in each region

% Check for NaNs in selected series
if any(any(isnan(Y))),
	error('NaN found in regional climate matrix Y')
end


temp1=ns2'; % convert to row vector -- number of stations each region
denom=temp1(ones(mc,1),:); % expand to matrix
Y = Y  ./ denom;  % climate variable averaged over specified
		% stations in region;  Y holds the regional-average climate
		% series to be used (indirectly) as the predictand variables in 
		% spatial recontruction. Note that Y(L1C,:) will be the calib-pd
		% rows of Y
[my,q]=size(Y);  % spatial model will have q predictands

%********************  STEP-3:SINGLE-SITE RECONSTRUCTIONS  *************
%
% 3.1 Check for validity of time periods 
% 3.2 Form predictor matrix
% 3.3 Estimate parameters of models
% 3.4 Reconstruct

% Form predictor matrix;  TL will have lag-0 data in cols

[TL,yrTL] = lagyr3(T,YRS(1,:),lg);
% Note that TL is identical to T.  This step was retained to
% keep parallel structure in spatrec1.m and spatrec0.m


TC = TL(L1T,:);  % subset of TL to be used for full calibration
[mTC,nTC]=size(TC);
% Initialize matrix to hold R-sq, adj R-sq and F-ratio for final eqn
a=NaN;
Cwgt=zeros(2,nt); % will hold regression coefs, constant first
STATS = a(ones(nt,1),ones(3,1));
DHS = a(ones(mTC,1),ones(nt,1)); % calib-pd recons
DHL = a(ones(n,1),ones(nt,1)); % entire-rows-of-T recons
II1 = a(ones(nt,1),:);
II2 = a(ones(nt,1),:);

%
% Loop over each tree
for i = 1:nt;
	c1 = zeros(2,1); % To hold regr coefficients for this tree
	I1 = [i ]; % for lag 0
	II1(i,:) = I1;
	y = D(L1C,i); % predictand
	if any(any(isnan(TC(:,I1))))
		error('A calib-period predictor value is NaN')
	end
	[I2,I4,stats,c,e,yhs]=stepr2(TC,y,I1);
	II2(i,1:length(I2))=I2;
	lenc = length(c);
	c1(1)=c(1);
	c1(I4+1)=c(2:lenc);
	Cwgt(:,i)=c1;
	STATS(i,:) = stats;
	DHS(:,i) = yhs;
	Tall = [ones(n,1) TL(:,I1)];
	DHL(:,i) = Tall * c1; 
end


% DHS now holds the single-site recons for the calib period
% DHL now holds recons for all years covered by T


%*********  SPATIALLY AVERAGE SINGLE-SITE RECONS *************

disp('Beginning To Spatially Average Single-Site Reconstructions')

L3=~isnan(IDX);
if n4 == 1; % special case: IDX Is not a matrix but a cv
	ns3=L3;
else
	ns3=   (sum((L3)'))'; % col vector of number of sites to average
end
if any(ns3<1), error('Some row of IDX is all NaN'), end

% Form logical matrix corresponding to IDX that can be used in 
% matrix multiplication 
L4 = ind2log(IDX,nt);
LX8 = isnan(DHS);
LX9 = isnan(DHL);
sum8 = sum(sum(LX8));
sum9=sum(sum(LX9));
DHSz=DHS;
DHLz = DHL;
DHSz(LX8)=zeros(sum8,1);
DHLz(LX9) = zeros(sum9,1);


X = DHSz * L4';  % tsm of sum of single-tree recons over sites
Xbig = DHLz * L4';

a = NaN;
X(LX8) = a(ones(sum8,1),:);
Xbig(LX9) = a(ones(sum9,1),:);

% Check for NaNs in selected series
if any(any(isnan(X))),
	error('NaN found in predictor matrix X')
end

temp1=ns3'; % convert to row vector -- number of sites in each group
denom = temp1(ones(mTC,1),:); % expand to matrix
X = X ./ denom;

denom=temp1(ones(n,1),:); % expand to matrix
Xbig = Xbig  ./ denom;


%************************ SPATIAL RECONSTRUCTIONS *****************
%
%******* Get col-index matrices for distinct recon-period data sets

YRS1 = [YRS(1,:); YRS(4,:); YRS(3,:)];
[L6,modnum]=uniqmods(Xbig,YRS1);
% L6 (mL6 x nL6)L pointer matrix telling which cols of X and Xbig
%		are active for each of mL6 models
% modnum (matrix)i which model (row of L6) applies for each
%		of the reconstruction years in Z as pointed to by 
%		row 2 of yrs
[nmods,dum1]=size(L6);  % nmods is number of different models

% Loop over each of the recon-year models, estimating parameters,
% getting reconstructed values, and saving model statistics

% Allocate 
anan=NaN;
mrec = sum(L2T); % number of reconstruction-period years
YP = a(ones(mrec,1),ones(q,1)); % reconstructed climate
EV=anan(ones(nmods,1),:);
RSQ1=anan(ones(nmods,1),ones(q,1));
RSQ2=anan(ones(nmods,1),ones(q,1));
NV=anan(ones(nmods,1),ones(4,1));

for n5=1:nmods;% 2;
	disp(['Beginning spatial model ',int2str(n5)])
	jcols =L6(n5,:);
	% Find out rows (years) in T to be reconstructed with this model
	L7 = modnum==n5;
	XXbig = Xbig(L2T,:);
	Xrec = XXbig(L7,jcols);
	Xcal = X(:,jcols);
	Ycal = Y(L1C,:);
	[yp,yrec,rsq1,rsq2,ev,nv]=osr(Xcal,Ycal,Xrec,mineigs,tcrit,kopt);
	NV(n5,:)=nv;
	RSQ1(n5,1:nv(2))= rsq1;
	RSQ2(n5,:)=rsq2;
	EV(n5) = ev;
	YP(L7,:) = yrec; 
end


end; % of kdope==1


disp('Beginning Crossvalidation')



%******************** CROSSVALIDATION -- SINGLE TREE  **********
%
% Must repeat the estimation of single-tree models using the same
% models as in the previous single-site reconstruction, but
% leaving out nout years from the estimation. the nout years are
% successively shifted by one and the model recalibrated
%
% Have nt different trees for which single-tree regressions are done.
% Regression models were previously estimated for these trees. We will
% use the same designated predictor variables as found in the full
% calibration.
%
% Outer loop over the the mTC leave-nout-out calibration periods.
% Inner over the nmods spatial models.
% 
% Within the outer loop, must re-calibrate the single-tree-models.
%
% Critical Critical output info includes, for each of the nummods models
% - mean, median, lo, hi R-squared for the leave-n-out calibrations
%		in single tree mods
% - r-squared and RE for the single-tree validations
% - mean EV for the leave-n-out spatial regressions
% - validation estimates of climate at the q points
% - RV computed from validation-period residuals
% 
%
% Keep these settings from the earlier steps
% - same predictor variables for the single-tree regressions
% - same number of eigenvectors retained on the predictands and
%    predictors (qq and pp from pcareg1.m).  Note that the
%    exact predictors of the pp are not set in stone. This seems
%    reasonable, because the same eigenvectors might not result
%    from PCA on the full calib period and on the full period
%    minus the years needed for the validation estimate
%
%  Currently available variables from earlier steps:
%
% TC (mTC x nTC)r lagged tree-ring data, mTC calib years,
%		nt*3=nTC cols (lags 0,-1,+1)
% II2 (nt x 3)i  cols of TC for models, NaN right-filled
% NV (nmods x 4)i [q qq p pp] number of predictands, retained
%	predictands, predictors, and retained predictors1 x nummods)i 

% Initialize key matrices
R2 = a(ones(nt,1),ones(mTC,1)); % Calib R-squared for the mTC
	% leave-3-out single-site models at each of the nt trees
EVcv=a(ones(mTC,1),ones(nmods,1)); % calibration stat
DHScv=a(ones(mTC,1),ones(nt,1)); % single-tree predictions of climate
Z = a(ones(nmods*mTC,1),ones(q,1)); % predicted climate at q points 
zgo = 1:mTC:(nmods*mTC-(mTC-1)); % starting row in Z of each 
	% mTC-row predicted submatrix
Mcv=Z; % calibration-period means needed for RE statistic


% Get the index for crossvalidation years
LTC=crospul2(mTC,0,0);



for k1 = 1:mTC;  % Loop over crossvalidation periods
	disp(['Starting crossv period ',int2str(k1)]);
	Lr = LTC(:,k1); % rows of calibration period to calibrate model on 
	yr5 = yr2(Lr); % year-cv for current leave-3-out calib period
	nobs = sum(Lr); % number of calibration years
	TC1 = TC(Lr,:); % predictor calibration data (lagged tree matrix)
	DHkey=a(ones(mTC,1),ones(nt,1));
	DHScv=a(ones(nobs,1),ones(nt,1));

	% Single-site regression, looping over sites
	for k2 = 1:nt;
		y = D(L1C,k2); % full calib-period predictand vector
		y = y(Lr);  % just the sub-period for crossvalid. fitting
		ii2 = II2(k2,:);
		npreds = sum(~isnan(ii2));
		ii2 = ii2(1:npreds);
		[bwt,r2,yhat] = regress1(y,[ones(nobs,1) TC1(:,ii2)]);
		 
		Tkey = [1   TC(k1,ii2)];
		DHkey(k1,k2) = Tkey * bwt; % prediction for the indep data
		R2(k2,k1)=r2; 
		%Tip: plot(R2); shows variability of R-squared from one
		% model to another as different groups of obs are omitted
		% [STATS(:,1) (R2(:,k))] shows comparison of calib R-squared for
		%		full-period calibation with calib-R-squared for
		%		kth leave-3-out calibation
		DHScv(:,k2) = yhat; % calib-period predictions for currrent
			% leave-3-out model
		% Tip: For time series plots of (1) observed single-site
		% climate series, (2) full-period calibration-period 
		% predictions, and (3) leave-3-out calibration-period
		% predictions for current model, do this:
		%  plot(yr2,D(:,3),yr2,DHS(:,3),yr5,DHScv(:,3))
	 end

	 
	% Spatially average single-tree reconstructions
	% Form logical matrix corresponding to IDX that can be used in 
	% matrix multiplication -- L4 formed previously
	LX10 = isnan(DHScv);
	dhkey = DHkey(k1,:);
	LX11 = isnan(dhkey);
	sum10 = sum(sum(LX10));
	sum11 = sum(LX11);
	DHScvz=DHScv;
	dhkeyz=dhkey;

	
	if sum10~=0,
		 DHScvz(LX10)=zeros(sum10,1);
	end
	if sum11~=0;
		 dhkeyz(LX11)=zeros(sum11,1);
	end

	% Weight single-site predictions over sites
	Xcv = DHScvz * L4';  % tsm of sum of single-tree recons over sites
	Xkey = dhkeyz * L4';

	a = NaN;
	if sum10~=0
		Xcv(LX10) = a(ones(sum10,1),:);
	end
	if sum11~=0,
		Xkey(LX11) = a(ones(sum11,1),:);
	end

	temp2=ns3'; % convert to row vector -- number of sites in each group
	denom = temp2(ones(nobs,1),:); % expand to matrix
	Xcv = Xcv ./ denom;
	Xkey = Xkey ./ temp2;

	% Check for NaNs in selected series
	if any(any(isnan(Xkey))),
		error('NaN found in predictor matrix Xkey')
	end

	if any(any(isnan(Xcv))),
		error('NaN found in predictor matrix Xcv')
	end




	% PCA REGRESSION FOR CROSSVALIDATION	
	

	for n5=1:nmods; % Loop over the different tree-group models
		disp(['   Starting tree-group model ',int2str(n5)]);
		jcols =L6(n5,:);
		Tcal = Xcv(:,jcols);
		Ycal = Y(L1C,:);
		Ycal = Ycal(Lr,:);
		Trec = Xkey(jcols);
		[ypcv,z1,rsq1,rsq2,ev,nv]=osr(Tcal,Ycal,Trec,mineigs,tcrit,kopt);
		%RSQ1cv(n5,1:nv(2))= rsq1;
		%RSQ2cv(n5,:)=rsq2;
		EVcv(k1,n5) = ev;

		 iz = zgo(n5) + k1-1; % row index in Z to store this prediction
		Z(iz,:) = z1; 
		Mcv(iz,:) = mean(Ycal); % calibration-period means needed for RE
  end



end;  % k1 loop over crossvalidation periods


% To compare time series of 
% the observed weighted climate variable at each tree-site
% with the single-site prediction based on independent data
%[(DHkey(k,:)  D(k,:)]
% plot(yr2,DHkey(k,:),yr2,D(k,:))





% Computer verification squared correlation coefficient and RE


RE = a(ones(nmods,1),ones(q,1)); % redcution of error at each gridpoint
RSQ3=RE;  % squared correl coeff between predicted and actual
	% climate for each spatial-model/gridpont
MSE=a(:,ones(q,1)); % initialize
MSM=MSE;

for i= 1:nmods ; % loop over spatial models
	Y1 = Z(zgo(i):(zgo(i)+mTC-1),:); % predicted climate
	M1 = Mcv(zgo(i):(zgo(i)+mTC-1),:)
	EM = Y - M1; % actual minus calib-period mean
	EY = Y - Y1; % residuals = observed minus predicted
	MSE = mean (EY .* EY); % rv of mean square errors at the 
	MSM = mean(EM .* EM);

	% RE statistic
	RE(i,:) = 1 - (MSE ./ MSM);

 	% Correlation coefficients
	for j = 1:q;
		yreal=Y(:,j);
		y1 = Y1(:,j);
		rr = corrcoef([y1 yreal]);
		RSQ3(i,j) = rr(1,2) .^2; 
	end
	
end

% note:
% For RE, will need calib-period means of predictand for each
% leave-3-out calib period and each spatial model
% This could be ((mTC x nmods) x q).  Then incomputing RE, would
% need to who knows what.
