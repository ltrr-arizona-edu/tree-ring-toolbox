function [YP,Z,RSQ1,RSQ2,EV,nv,hf]=osr2(X,y,X1,mineigs,tcrit,kopt)
% Use: [YP,Z,RSQ1,RSQ2,EV,nv]=osr2(X,y,X1,tcrit,kopt)
% osr2.m  -- orthogonal spatial regression, without forcing to z-scores (as in
%    osr.m, and handling only single predictand
%
% D Meko 10-18-97
%
%
%*************** IN ARGS **********************************
%
% X (mm x p)r  tsm of predictors. Typically scaled so that
%		differences in variance important, as in series
%		produced by reglag1.mstructions
% y (mm x 1)r cv time series of predictand
% X1 (mX1 x p)r set of predictors for reconstruction years
% mineigs (1 x 1)i min number of predictor variables for
%		which eigenvalue discarding is used
% tcrit (1 x 1)r significance level (typically 0.95) for
%	  including predictors in pcareg1.m
% kopt(1 x 2)i options
%		kopt(1)-- selection of important predictor eigs
%			1 -- use all
%			2 -- eigenvalue of 1
%		kopt(2) -- reserved
%		kopt(3) -- predictor selection option
%			1 -- use all predictors
%			2 -- only those whose reg coeffs sig at tcrit%
%******************  OUT ARGS ******************************
%
% YP (mm x 1)r reconstructed predictand for calib period
% Z (mX1 x 1)r reconstructed predictand for data X1
% RSQ1 (1 x ?)r R-squared for each non-discarded predictand eigenv
% RSQ2 (1 x q)r R-squared for each point or station in y
% EV (1 x 1)r summary R-squared computed from RSQ1 and pctg of
%		predictand variance accounted for by each predictanc eign
% nv(1 x 5)i number of variables:
%	 1 - q=full number of predictand eigenvectors
%  2 - qq=number of those used (after discarding unimportant eigs)
%  3 - p=full number of predictor eigenvectors
%  4 - pp=number of potential predictors after discarding eigenvectors
%  5 - number of predictors actually used after t-test deletion
% hf(mX1 x q)L  extrapolation (1) vs interpolation(0) flag for
%		reconstructed data Z
%
%****************** USER-WRITTEN FUNCTIONS CALLED ***********
%
% eigen1.m -- PCA on predictors and predictands
% mce1.m -- minimum covering ellipsoid to flag extrapolations
% pcareg1.m
%
%********************  NOTES ******************************
%
% Method follows Cook, Edward R., Briffa, Keith R., and 
% Jones, Philip D., 1994.  Spatial regression methods in 
% dendroclimatology:  a review and comparison of two techniques.
% International Journal of Climatology 14, 379-402.
%
% t-test for significance of PC-regression coefficients: 
% Mardia, K. V., Kent, J. T, and Bibby, J. M., 1979.  Multivariate
% Analysis.  Academic Press Limited, 518 pp.  See p. 244-246.
%
% Interpolation vs extrapolation:
% Weisberg, Sanford, 1985.  Applied Linear Regression.  
% John Wiley & Sons, 324 pp.  See p. 235-237.


%
% STEPS
% - convert y to z-scores, using calib period mean and std deviation
% - PCA on covariance mtx of X  and delete unimportant components to reduce
% - Convert retained PC scores to z-scores
% - Regress predictand on predictor amplitudes, further
%	 simplifying by discarding predictors that fail t-test for 
%   significance
% - Transform reconstruction equation to give predictions of original
%	 climate variables
% - Diagnostic statistics


%******************  Initialize, size, preallocate *************

[mm,p]=size(X);
[dum1,q]=size(y);
if mm~=dum1
	error('X and y must be same row-size');
end
if q~=1;
	error('y must be column vector');
end


[mX1,nX1]=size(X1);
if nX1~=p,
	error('X1 must have save col size as X')
end

%************  Compute means and standard deviations
Xmean=mean(X);
ymean=mean(y);
Xstd = std(X);
ystd = std(y);

% Compute z-score time series
Xz=zscores(X);
yz=zscores(y);

% Expand mean and std deviation of predictand into col vectors
yym = ymean(ones(mm,1),:);
yys = ystd(ones(mm,1),:);


%******************  PCA for data reduction ***************
%
% Note that the scores are compute by "F*V" within eigen1.m not "F*B"
% The scores will thus not be unit standard deviation
%
%
[RE,E,LE,SE,U,CE]=eigen1(X,2); % PCA on covariance matrix of X


if kopt(1)==1; % use all predictors
	pp=p;
elseif kopt(1)==2; % eigenvalue of 1
	if p>mineigs,
		pp =  max([sum(diag(LE)>1) mineigs]);
	else
		pp=p;
	end
end



U = U(:,1:pp);

LEsub = diag(LE);
LEsub = LEsub(1:pp);
LEpp = LE(1:pp,1:pp);
%&&&&&&&  

%*************** Prior to regression, normalize each amplitude by the
% square-root of its eigenvalue


LLE = 1 ./ sqrt((diag(LEpp))');
LLE = diag(LLE);

Uz = U * LLE;

%**************  Estimate regression coefficients

if kopt(3)==1; % skip the t-test, use all potential predictors

	A = inv(Uz'*Uz)*Uz'*Vz;

elseif kopt(3)==2; % sig test on coefficients -- use 95%
	[atilde,t] = pcareg1(LEsub,Uz,y,E(:,1:pp),tcrit);
	A(:,1) = atilde;
end; % of if kopt(3)

B = E(:,1:pp)* A * (F(:,1:qq))';
nused = sum(A~=0);

nv = [q qq p pp nused];

% ***** FRACTIONAL VAIANCE OF DERIVED qq PREDICTANDS EXPLAINED 
% Eqn 16 in Cook et al.; eqn 8.8.6 in Mardia et al.
RSQ1 = (diag(A' * A))'; % rv of R-squared for each score-predictand


%*************** FRACTIONAL VARIANCE OF ORIGINAL q PREDICTANDS
%  (Cook et al, top of p. 385)
EV = trace(A'*A*LF(1:qq,1:qq))/q; % average fractional explained variance
	% in terms of the original q predictands




%************  Calibration period predicted values of amplitudes

VPz = Uz*A;
VP = VPz * inv(LLF);

% Cal-period predictions expressed as original variables
% Note that eigen1.m converts the data to z-scores before computing
% eigenvectors.  We must consider this in the back-transformation.
% First get the prediction of standardized departures of y in terms
% of the 


yPz = VP * (F(:,1:qq))';

% Now get the prediction in terms of a coefficient matrix
% to postmultiply the "standardized" original X variables

B =E(:,1:pp)*LLE *A * inv(LLF) * (F(:,1:qq))';
%yPz2 = Xz * B;



%****** SPATIAL DISTRIBUTION OF R-SQUARED  ************************
%
err = yz - yPz;
RSQ2 = 1.0 - (var(err) ./ var(yz));



YP = YPz .*YYs + YYm; % predicted Y


%********************* Reconstruction-period predicted amplitudes

X1m = Xmean(ones(mX1,1),:); % expand calib-period means to size of X1
X1s = Xstd(ones(mX1,1),:); % expand calib-pd std devs to size of X1
XX1 = (X1-X1m) ./ X1s; % XX1 standardized by calib-period mean, std dev

Z1 = XX1*B;
ZZm = Ymean(ones(mX1,1),:); 
ZZs = Ystd(ones(mX1,1),:);

Z = Z1 .* ZZs + ZZm;


%*****  Interpolation vs Extrapolation *************************
%
% Calibration-period predictor data for the regression is Uz, 
% which is zero-mean. Let the corresponding transformed data for the
% years outside the calibration period be U1z:
%
% Get used subset of columns of E
% Find out which of the pp predictor amplitudes were used in any of the
% regressions to predict the qq predictand amplitudes
j1 = A~=0;  % logical col vector to used cols of E(:,1:pp)
E1 = E(:,1:pp);
E1 = E1(:,j1);

U1z = XX1*E1;

% Get matrix of used calibration-period predictor eig-amps
U2z = Uz(:,j1);

% Matrix of corresponding data outside the calibration period
% is already U1z 

% Stack non-calibration and calibration predictor data 
U3z=[U1z;U2z];

% The call for minimum covering ellipsoid is
% [S,hmax]=mce1(X,yrs1,yrs2), where X is U3z, yrs1 are start and end
% year for X, and yrs2 are for calibration period
yrs1=[1 (mX1+mm)];
yrs2=[(mX1+1) (mX1+mm)];
[S,hmax]= mce1(U3z,yrs1,yrs2);
hf=S(1:mX1,3); % flag for outside-calibration-period years
	% 1 = extrapolation, 0 =interpolation
