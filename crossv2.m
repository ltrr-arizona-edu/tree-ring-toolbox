function [V1,V2]=crossv2(X,y,nkeep)
% crossv2: exploratory cross-validation of a multiple recflow1.m reconstruct models
% CALL: [V1,V2]=crossv2(X,y,nkeep);
%
% Meko 10-15-97
%
%****************** IN
%
% X (mX x nX); % predictor data (for calibration period)
% y (my x 1); % predictand data for same years
% nkeep (1 x 1)i  % max  number of PCs to keep as potential predictors
%
%
%*************** OUT 
%
% V1 (nkeep x nkeep)r cross-validation rmse for various combinations of
%		number of PCs retained as potential predictors (cols), and number
%		actually usd as predictors
%
% V2 (nkeep x nkeep)r likewise, but reduction of error statistic
%
%************** NOTES
%
% Method.  Fits successively more complicated PCR models.  Number of PCs
% retained is varied from 1 to nkeep.  Within that outer loop, number of
% PCs used for prediction is varied from 1 to the maximum possible.  
% Predictors are added in  order of the size of their standardized regression
% coef.


a=NaN;
pp=nkeep; % rename for convenience in coding

% Get length of calibration period, and number of predictor variables
[nyrs1,nvar]=size(X);

% Check that y correct length
if size(y,1) ~= nyrs1;
	error('X and y must be same row size');
end


%************  Compute full-calibration period means 
Xmean=mean(X);
ymean=mean(y);


% Expand means into vectors or matrices
yym = ymean(ones(nyrs1,1),:);
XMN = repmat(Xmean,nyrs1,1);


%****************  CROSS-VALIDATION BY LEAVE-1-OUT METHOD


% Allocate
ykey = a(ones(nyrs1,1),ones(nkeep,1),ones(nkeep,1)); % cross-validation estimated time series 
V1 = a(ones(nkeep,1),ones(nkeep,1));  % rmse for cross-valid
V2 = V1; % reduction of error

C1 = ykey; % R-squared for calibration
C2 = ykey; % root-mean-square error for calibration
C3 = ykey; % re statistic for calibration


%-------- Logical matrix of observations to use for each cross-v model
IX=crospul2(nyrs1,0,0); % the zeros mean no neg lags or pos lags in model
IX=logical(IX);
% each col of IX goes with a cross-validation model
% rows of IX correspond to years of y and X
% entry of IX is 1 if observation to be used in fitting model, 0 otherwise

mm=nyrs1-1; % Number of years to calibrate on

% Loop over the calibration period years
for n = 1:nyrs1;
	ix = IX(:,n); % row pointers for this model
	% Get subset of X and y data for calibration
	XX = X(ix,:);
	yy = y(ix);

	% Get the single row of X data to be used for cross-valid estimate 
	Lzero=~ix;
	if sum(Lzero)~=1; 
		error('Should have only one zero element in ix here');
	end
	xkey = X(Lzero,:);

	% Compute calibration-period means; expand mean; mean-correct data
	XXmean=mean(XX);
	yymean=mean(yy);
	yyym = yymean(ones(mm,1),:);
	XXMN = repmat(XXmean,mm,1);
	XXcent =XX - XXMN;
	yycent = yy -yyym;

	% Also center the row of predictor data to be used for cross-val estimate
	xkeycent=xkey-XXmean;
		
	%-----------PCA on covariance matrix 
	[RE,E,LE,SE,U,CE]=eigen1(XXcent,2); 

	for k = 1:nkeep; % Loop over number of PC's retained as potential predictors

		%------Reduce the matrices of PC scores and eigenvalues as specified by k
		Usub = U(:,1:k); % PC scores, reduced set
		LEsub = diag(LE); % eigenvalues as a col vector
		LEsub = LEsub(1:k); % the reduced set of eigenvalues as a col vector
		LEpp = LE(1:k,1:k); % the reduced set of eigenvalues as a matrix

		%-----------  Estimate regression coefficients on retained PC scores
		atilde = inv(Usub'*Usub)*Usub'*yycent; % with all predictors in model

      %--------- Sort regression coeffs from least significant to most
      asize = atilde .* sqrt(mm*LEsub);
		[s,js]=sort(abs(asize)); % sorted smallest to largest

		for j = 1:k; % Loop over number of actual predictors used in model

			% Must change all except the j most significant atilde to zero
			asub=atilde; % copy atilde because must change within this inner loop
			nzero = length(atilde)-j; % change this many to zero
			jzero = js(1:nzero); % these are the elements of atilde to change to zero
			asub(jzero)=0;

			%------------ Express as coefficients on centered original predictors
			B = E(:,1:k)* asub; % regression coefs on the original predictor variables

			%-------Prop var explained ( see ewn 8.8.6 in Mardia et al. 1979)
			rsqeach = (mm * LEsub .* (asub .* asub)) / (yycent' * yycent) ; % cv of
					%R-squared each component
			C1(n)=sum(rsqeach);  % total proportion of variance explained in calibration

			%------------ Calibration period predicted values
			yyhat = Usub*asub + yyym;

			%------------ Regression residuals
			ehat=yy-yyhat; 

			%----------------- Calibration rmse and RE (bogus RE)
			[mae,rmse,re]=rederr(yymean,yymean,yyhat,yy);
			C2(n)=rmse(1); % this is root-mean-square error of reconstruction for calib
			C3(n) = re(1); % reduction of error statistic for calibration; based on calib
				%		period mean of predictand as the null reconstruction

			% Estimate the predictand for the left-out observation
			ykey(n,k,j) = (xkeycent * B) + yymean; 
		end
	end
end


%----------- CROSS VALIDATION  STATISTICS

for k = 1:nkeep; % loop over number of PCs retained
	for j = 1:k; % and over number used in model
		ygen = ykey(:,k,j); % the series generated by cross-validation
		[mae,rmse,re]=rederr(ymean,ymean,ygen,y);
		V1(j,k) = rmse(1);  % root-mean-square error of cross-validation estimates
		V2(j,k) = re(1); % reduction of error stat, based on comparison of recon with
			%		a null reconstruction using the calibration period mean as the 
			%		value in each year
	end
end
V1=flipud(V1);
V2=flipud(V2);
