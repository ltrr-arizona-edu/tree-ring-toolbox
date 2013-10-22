function [v1,v2,v3]=crossv1(X,y,nkeep,npred)
% crossv1: cross-validation of a single recflow1.m reconstruction model
% CALL: [v1,v2]=crossv1(X,y,nkeep,npred);
%
% Meko 10-15-97
%
%****************** IN
%
% X (mX x nX); % predictor data (for calibration period)
% y (my x 1); % predictand data for same years
% nkeep (1 x 1)i  % number of PCs to keep as potential predictors
% npred (1 x 1)i  % number of potential predictors to use in reconstruction
%		model (see notes)
%
%
%*************** OUT
%
% v3 (my x 1)r   cross-validation residuals
%
%************** NOTES
%
% npred.  the npred potential predictors are those with the highest
%		standardized regression coefficients

a=NaN;
pp=nkeep; % rename for convenience in coding

% Get length of calibration period, and number of predictor variables
[nyrs1,nvar]=size(X);
mm=nyrs1;

% Check that y correct length
if size(y,1) ~= nyrs1;
	error('X and y must be same row size');
end


%************  Compute full-calibration period means 
Xmean=mean(X);
ymean=mean(y);


% Expand means into vectors or matrices
yym = ymean(ones(mm,1),:);
XMN = repmat(Xmean,mm,1);


%****  'center' predictand and predictor data  (see Marsdia, p. 244, eqn 8.8.1
Xcent = X - XMN;
ycent = y -yym;


%****************  CALIBRATE MODEL USING ALL CALIBRATION YEARS

%-----------PCA on covariance matrix of centered X
[RE,E,LE,SE,U,CE]=eigen1(Xcent,2); 


%------  Reduce the matrices of PC scores and eigenvalues as specified by nkeep
U = U(:,1:pp); % PC scores, reduced set
LEsub = diag(LE); % eigenvalues as a col vector
LEsub = LEsub(1:pp); % the reduced set of eigenvalues as a col vector
LEpp = LE(1:pp,1:pp); % the reduced set of eigenvalues as a matrix

%-----------  Estimate regression coefficients on retained PC scores
atilde = inv(U'*U)*U'*ycent; % with all predictors in model

% Must change all except npred most significant elements of atilde to zero
% Significance is proportional to value of atilde times the 
% sqrt (sample size * eigenvalue) (see Mardia,  p. 245, eqn 8.8.5)
asize = atilde .* sqrt(nyrs1*LEsub);
[s,js]=sort(abs(asize)); % sorted smallest to largest
nzero = length(atilde)-npred; % change this many to zero
jzero = js(1:nzero);
atilde(jzero)=0;

%------------ Express as coefficients on centered original predictors
B = E(:,1:pp)* atilde; % regression coefs on the original predictor variables

%-------Proportion of variation expd ( see ewn 8.8.6 in Mardia et al. 1979)
rsqeach = (mm * LEsub .* (atilde .* atilde)) / (ycent' * ycent) ; % cv of
			%R-squared each component
rsq=sum(rsqeach);  % total proportion of variance explained

%------------ Calibration period predicted values
yhat = U*atilde + yym;

%------------ Regression residuals
ehat=y-yhat; 

%----------------- Calibration rmse and RE (bogus RE)
[mae,rmse,re]=rederr(ymean,ymean,yhat,y);
rmse1=rmse(1); % this is root-mean-square error of recons for full model
re1 = re(1); % reduction of error statistic for calibration; based on calib
%		period mean of predictand as the null reconstruction





%****************  CROSS-VALIDATION BY LEAVE-1-OUT METHOD


% Allocate
ykey = a(ones(nyrs1,1),:); % cross-validation estimated time series 
C1 = ykey; % R-squared for calibration
C2 = ykey; % root-mean-square error for calibration
C3 = ykey; % re statistic for calibration




%-------- Logical matrix of observations to use for each cross-v model
IX=crospul2(nyrs1,0,0); % the zeros mean no neg lags or pos lags in model
IX=logical(IX);
% each col of IX goes with a cross-validation model
% rows of IX correspond to years of y and X
% entry of IX is 1 if observation to be used in fitting model, 0 otherwise


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



	mm=nyrs1-1; % number of years to calibrate on

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

	%------Reduce the matrices of PC scores and eigenvalues as specified by nkeep
	U = U(:,1:pp); % PC scores, reduced set
	LEsub = diag(LE); % eigenvalues as a col vector
	LEsub = LEsub(1:pp); % the reduced set of eigenvalues as a col vector
	LEpp = LE(1:pp,1:pp); % the reduced set of eigenvalues as a matrix

	%-----------  Estimate regression coefficients on retained PC scores
	atilde = inv(U'*U)*U'*yycent; % with all predictors in model
   % Must change all except most significant npred of atilde to zero
   asize = atilde .* sqrt(mm*LEsub);
	[s,js]=sort(abs(asize)); % sorted smallest to largest
	nzero = length(atilde)-npred; % change this many to zero
	jzero = js(1:nzero);
	atilde(jzero)=0;

	%------------ Express as coefficients on centered original predictors
	B = E(:,1:pp)* atilde; % regression coefs on the original predictor variables

	%-------Proportion of variation explained ( see ewn 8.8.6 in Mardia et al. 1979)
	rsqeach = (mm * LEsub .* (atilde .* atilde)) / (yycent' * yycent) ; % cv of
		%R-squared each component
	C1(n)=sum(rsqeach);  % total proportion of variance explained in calibration

	%------------ Calibration period predicted values
	yyhat = U*atilde + yyym;

	%------------ Regression residuals
	ehat=yy-yyhat; 

	%----------------- Calibration rmse and RE (bogus RE)
	[mae,rmse,re]=rederr(yymean,yymean,yyhat,yy);
	C2(n)=rmse(1); % this is root-mean-square error of reconstruction for calib
	C3(n) = re(1); % reduction of error statistic for calibration; based on calib
	%		period mean of predictand as the null reconstruction

	% Estimate the predictand for the left-out observation
	ykey(n) = (xkeycent * B) + yymean; 
end


%----------  Cross-validation statistics
[mae,rmse,re]=rederr(ymean,ymean,ykey,y);
v1 = rmse(1);  % root-mean-square error of cross-validation estimates
v2 = re(1); % reduction of error stat, based on comparison of recon with
%		a null reconstruction using the calibration period mean as the 
%		value in each year
v3=y-ykey;  % cross-validation residuals
