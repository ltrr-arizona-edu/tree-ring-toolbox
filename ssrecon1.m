function [yhat,b,S,SC,SCV]=ssrecon1(YB,XB,nX,IY,IX,yrs,lg,maxlag,lagmask,cvopt)

% Single-site reconstruction at a network of gridpoints.

% D. Meko 11-20-92
%  Revised 3-24-94 to allow masking of lags by lagmask arg

%***********   USER WRITTEN FUNCTIONS NEEDED --- NONE
%
% crospull.m  set up index array for crossvalidation
% sos.m  multiple linear regression and sum of squares calcs

%***********   INPUT ARGS
%
% YB (m1 x n1) predictand array, mYB years, nYB variables
% XB (m2 x n2) lagged-predictor array, assumed to be generated previously
%		by lagyr3.m
% nX (1x1) number of unlagged predictor sites (some fraction of ncols
%		in IX)
% IY (m3 x 1) index to cols of YB specifying variables to reconstruct
% IX (m4 x 1) index to cols of XB specifying matching predictor variables.
%    Note that m3=m4, and that m3 can be smaller than n1.
% yrs (4x2) beginning and ending years of :
%		row 1: YB
%		row 2: XB
%		row 3: calibration period
%		row 4: entire reconstruction
% lg  (1x2) number of neg and pos lags used in building XB with lagyr3.m
% maxlag (1 x 2) maximum number of -,+ lags to be used in the regression
% lagmask(1 x (1 + sum(lg))   mask allowing some lags in lagged predictor
%		matrix to be excluded as predictors in model.  By example, if
%		lg=[1 1] you have specified that XB has lags 0,-1,1, reading
%		sets of columns left to right. 
%       If no + or - lags are to be included as predictors,
%		A compatible setting is lagmask=[1 0 0], which effectively
%		deletes from consideration the -1 and +1 lags.  The function will
%		bomb if you have not thought this through.  lagmask contains only
%		as many "1" values as predictors you want to use in the model.
%	cvopt (1x1) 'Y'=do cross validation,   'N'=no
%************  OUTPUT ARGS
%
% yhat (m5 x n5) the long-term reconstructed series.  Note that
%		n5=m3, and m5 is possibly as large as m2 if no lags or delays.
% b (m3 x (2 + sum(lags)) the regression coefficients for each model
%		First the constant term, then the zero lag, then lags -1, -2 ...
%		then lags +1, +2,...
% S (m3x3) R-squared, adjusted R-squared, and F-level for regr eqn
% SC (m3x6) summary statistics for calibration data (full model)
% SCV (m3x6) corresp statistics for crossvalidation. Cols defined as:
%		1. r-squared, squared Pearson r between observed and predicted
%		2. mean absolute error using reconstruction as estimates
%		3. mean absolute error using calib-period observed mean as estimates
%		4. rmse corresponding to (2)
%		5. rmse corresponding to (2)
%		6. reduction of error statistic


%*************  YEARS AND ROWS CALCULATIONS

yr1= (yrs(1,1):yrs(1,2))';  % years vector for YB
yr2= (yrs(2,1):yrs(2,2))';  % years vector for XB
yr3= (yrs(3,1):yrs(3,2))';  % years vector for calib period
yr4= (yrs(4,1):yrs(4,2))';  % years vector for entire reconstruction

% Make logical pointers to help in pulling out years

L1YB= yr1>=yrs(3,1) & yr1<=yrs(3,2);  % to calib pd years in YB
L1XB= yr2>=yrs(3,1) & yr2<=yrs(3,2);  % to calib pd years in XB
L2XB= yr2>=yrs(4,1) & yr2<=yrs(4,2);  % to long-term recon pd in XB
L1yh= yr4>=yrs(3,1) & yr4<=yrs(3,2);  % to calib-pd years in yhat

%*************  PREALLOCATE, SIZE, INITIALIZE

[m1,n1]=size(YB);
[m2,n2]=size(XB);
m3=length(IY);
m4=length(IX);


% check that lagmask proper size and only 1s or 0s
if length(lagmask) ~= (sum(lg)+1), error('lagmask wrong length'), end;
onelm=ones(1,length(lagmask));
zerlm=zeros(1,length(lagmask));
LLL=lagmask~=onelm;
MMM=lagmask~=zerlm;
if any(LLL & MMM), error('lagmask must contain zeros or ones'), end;

a=NaN;
yhat=a(ones(length(yr4),1),ones(length(IY),1));
yhat1=a(ones(length(yr3),1),ones(length(IY),1));
ncoefs=1+sum(lagmask);  % maximum possible number of coefs
b=zeros(m3,ncoefs);  % fill reg coef array with zeros
SC=a(ones(m3,1),ones(6,1));
SCV=a(ones(m3,1),ones(6,1));
ycross=a(ones(length(yr3),1),:);


Jsize=1+lg(1)+lg(2);
JB=a(ones(n1,1),ones(Jsize,1));

S=zeros(m3,3); % will hold R-squared, adj R-sq, and F-ratio for equations

%****************  BUILD AN INDEX ARRAY TO COLS OF XB

row1=[1:nX:(1+sum(lg)*nX)]; % First row of JB
JB=row1(ones(n1,1),:);  % dupe that row
cv1=(1:n1)'-1;  % col vector
arr1=cv1(:,ones(Jsize,1));  % dupe the col vector
JB=JB+arr1;  % This is the index array.  Row 1 of JB hold the
%	col subscripts of XB corresp to appropriate variables for 
%  model model for variable in col 1 of YB, etc.



%************* POSSIBLY DELETE ONE OR MORE COLS FROM JB TO RESTRICT LAG
% MODEL


% lagmask masks out selected lags.  In setting lagmask, you need to know
%  that col 1 of JB corresps to zero lag variable
%  col 2 ... to -1 lag
%  col 3 ... to -2 lag etc
%  col 2+lg(1) to  +1 lag etc
%  So for example, JB has 3 cols if lag -1 and +1
%  and passing lagmask=[1 1 0] will result in picking only the lag-0 and lag -1
%  predictors.  By default, fill lagmask with "1"s, which will omitt no cols.

JB=JB(:,lagmask);  % cull out the desired lags of predictors



%*****************    DO THE REGRESSIONS

for i = 1:m3;  % Loop for each predictand (each model)
	Y=YB(L1YB,IY(i));  % predictand cv
	i1=IX(i);
	jb=JB(i1,:);
	X=XB(L1XB,jb);  %  predictors 
	[bsub,Ssub]=sos(Y,X);  % call function to do regression and
	S(i,:)=Ssub;
	b(i,:)=bsub';
end


%*****************  RECONSTRUCT

% Have yr4 (cv) with years for reconstruction
% Have L2XB as logical pointer to those years in XB
cvones=ones(length(yr4),1);  % cv of ones, same length as reconstruction

for i=1:m3;  % Loop for each predictand
	i1=IX(i);  % which predictor variable for this predictand
	jb=JB(i1,:); % corresp col subscripts in XB for zero-lag and other
%		lags
	X=XB(L2XB,jb); % Predictor variables
	yhat(:,i)=[cvones X] * (b(i,:))';
end


%******************  CROSSVALIDATE

if cvopt=='Y'
	onesy=ones(length(yr3),1);  % ones vector for calib period
	yhat1=yhat(L1yh,:);  %  subset of reconstruction array covering calib pd
	ipull=crospull(length(yr3),lg(1),lg(2));  % logical cv whose cols
%		corresp to pointers to be used in crossvalid for a specific year. 
%	   Zeros indicate which observations to be ommitted.
	for i=1:m3;  % Loop for each predictand
		YY=YB(L1YB,IY(i));
		YYh= yhat1(:,i);  % reconstructed from full model
		i1=IX(i);
		jb=JB(i1,:);
		XX=XB(L1XB,jb);
		for j=1:length(yr3);  %loop for each year
			X=XX(ipull(:,j),:);
			Y=YY(ipull(:,j),:);
			
			[bsub,Ssub]=sos(Y,X);
			
			ycross(j)=[1 XX(j,:)] * bsub;  % estimate for the key leftout yr
		end

		rr1=(corrcoef(YYh,YY)) .^2;  % equiv to R-squared
	     rr2=(corrcoef(ycross,YY)) .^2; % sqd correl of cross-val est and obs
		mean1=mean(YY);  % observed mean for calib period
		[mae,rmse,re]=rederr(mean1,mean1,YYh,YY);  % full-model estimate vs
			% observed
		SC(i,[1 2 3 4 5 6])=[rr1(1,2) mae(1) mae(2) rmse(1) rmse(2) re(1)];

		[mae,rmse,re]=rederr(mean1,mean1,ycross,YY);  % cross-val estimate vs
			% observed
		SCV(i,[1 2 3 4 5 6])=[rr2(1,2) mae(1) mae(2) rmse(1) rmse(2) re(1)];

		
	end	


end
	
