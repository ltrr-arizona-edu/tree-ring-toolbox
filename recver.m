function [c,s1,s2,s3,yhat]=recver(b,Ib,yxf,yrs,k1);
% Reconstruction and verification

% D. Meko, 4-22-92

% Applies regression coefficients to a previously-built array of
% predictand and lagged predictors to generate reconstruction.  
% Evaluates fit of observed and model data in specified calibration and
% verification periods.  Period specified for verification should not 
% overlap with calibration period.

% Regression coefficients b, and their pointer vector Ib are assumed
% to have been estimated elsewhere (e.g., BMDP)

%******  INPUT ARGS
%
% b (n2 x 1)  Estimated regression coefficients, constant term first
% Ib (n3 x 1) Pointer to cols of yxf telling which vars the first
%     n3=n2-1 coefs are to be applied to.
% yxf (m1,n1) predictand/lagged-predictor array, formed by lagyr2.m
%		Col 1 is year, col 2 is predictand, next cols are unlagged 
%		predictors, next are negative lagged predictors, then positive
%		lagged predictors (see comments in lagyr2.m)
% yrs (2 x 2) begin and end years of  calib period (row 1) and 
%		verif period (row 2)
% k1 (1 x 1) option: k1=1,  full verification/reconstruction
%                    k1=2,  recon only, with R-squared stat  for check


%********  OUTPUT ARGS
%
% c (4 x 7)  comparison stats for actual and recon data
%  row 1: calib period, actual data
%  row 2: calib period, recon  data
%  row 3: verif period, actual data
%  row 4: verif period, recon  data
%
%       col 1 -- mean           col 5 -- max
%       col 2 -- median         col 6 -- r1 (first order autocorrel)
%       col 3 -- st dev         col 7 -- N, sample size
%       col 4 -- min
% 
% s1 (3 x 2) sign-test results for verification period (p. 330, Fritts1976)
%	row1: total number of cases N 
%	row2: minimum of number of agreements or disagreements
%	row3: significance level: .95 , .99, or 0
%
%		col 1 --  test on departures
%		col 2 --  test on first difference series
%
% s2 (2 x 9) reduction of error and related statistics
%		row 1:  computations for verification period
%		row 2:  computations for calibration period
%
%			col 1: sample size
%			col 2: MAE of reconstruction
%			col 3: MAE if calib-period mean were reconstruction
%			col 4: MAE if verif-period mean (obs data) were reconstr.
%			col 5: RMSE of reconstruction
%			col 6: RMSE if calib-pd mean were reconstruction
%			col 7: RMSE if verif-pd mean (obs data) were reconstruction
%			col 8: RE   if calib-pd mean were reconstruction
%			col 9: RE   if verif-pd mean (obs data) were reconstruction
%
%
% s3 (2 x 8) Correlation statistics  (zeros for  N/A)
%	row 1: calibration period
%	row 2: verification period
%
%		col 1:  Coef of determination (R) for calib;  correl coef (r) for 
%						verification period
%		col 2:  R-squared and r-squared
%		col 3:  first-order autocorrelation coefficient, r1
%		col 4:  F-statistic for R-squared, or
%				t-statistic for r
%		col 5:  95% prob point of F distribution
%     		95% prob point of t distribution
%		col 6:  99% prob point of F distrib
%				99% prob point of t distrib
%		col 7:  sample size, unadjusted
%		col 8:  sample size, adjusted  = N * (1-r1)/(1+r1)
%
% yhat (m1 x 2)  reconstruction, year as col 1


%*****  OTHER USER FUNCTIONS NEEDED
%
% signtest.m  --  does sign test
%	signtbl.mat -- automatically loaded table of sig points
% rederr.m  -- reduction of error and related statistics
% ftable95.mat -- 95% points for F test

%********  Check args;  size arrays;  initialize arrays


if ~exist('f95'), load ftable95, end;  % table for F-test
if nargin ~= 5, error('NEED 5 INPUT ARGS');  end;
if nargout ~= 5, error('NEED 3 OUTPUT ARGS');  end;

[m1,n1] = size(yxf);
n2=length(b);
n3=length(Ib);
[m4,n4]=size(yrs);

c=zeros(4,7);  s3=zeros(2,8); s2=zeros(2,9); yhat=zeros(m1,2);

yrs(3,1)=yxf(1,1); yrs(3,2)=yxf(m1,1);  % get start, end years of yxf
k4='N';  % initial setting for keyboard mode for plots

%******  Echo-check of input

clc; home;
disp('Echo chck of input args, including first 5 years, cols of yxf');
disp(' ');
disp('b transpose'); disp(b');  disp('Ib transpose');  disp(Ib');
disp('yrs'); disp(yrs);  disp('k1');  disp(k1);
disp(' ')
disp('Press any key for rest of echo-check')
pause
clc; home;
disp('FIRST 5 YEARS OF ARRAY yxf, TRANSPOSED')
disp('')
disp((yxf(1:5,1))')
disp((yxf(1:5,2))')
disp((yxf(1:5,3))')
disp((yxf(1:5,4))')
disp((yxf(1:5,5))')
disp('')
disp('Any key to continue');
pause
clc
%*******  Some consistency checks on input
% Ib should have one fewer element than b.
% Calib and verif periods must lie within years of yxf.

if n3~=n2-1, error('length(Ib) ~= length(b)-1'); end;

L1 = [yrs(1,1)<yrs(3,1)  yrs(2,1)<yrs(3,1)  yrs(1,2)>yrs(3,2) ...
yrs(2,2)> yrs(3,2)];
if any(L1),  error('yxf years do not include calib and verif pds'); end;

%*********  Form year vectors and see if calib and verif pds overlap
%   Overlap is not a problem if mode is reconstruct only (k1=2)

yr1=(yrs(1,1):yrs(1,2))';
yr2=(yrs(2,1):yrs(2,2))';
yr3=(yrs(3,1):yrs(3,2))';

L2c= yr3 >= yrs(1,1) & yr3 <= yrs(1,2);
L2v= yr3 >= yrs(2,1) & yr3 <= yrs(2,2);

L2s=sum([(L2c)';  (L2v)']);
if (max(L2s) >1) & k1 ==1, 
	error('CALIB AND VERIF PERIODS OVERLAP!');
end


% Compute yhat (the reconstruction), and a value for R-squared for the
% calibration period for double-checking against BMDP output.

X= [ones(m1,1)  yxf(:,Ib)];  % sub-matrix of selected predictors, and
%       a "ones" column for the constant
yhat = [yxf(:,1)  X*b];
rtemp=corrcoef(yhat(L2c,2),yxf(L2c,2));
s3(1,1)=rtemp(1,2);  % correl between actual and recon for calib period
clc; home;
disp(['Calibration period R-squared = ',num2str(s3(1,1) ^2)]);
disp('')
disp('Above value should match BMDP output.');
disp('')
disp('Press any key to continue');
pause


if k1 == 1 ;  %  want full verification info

%  Begin by filling the c array

	yc1= yxf(L2c,2);
	yc2= yhat(L2c,2);  % calibration period data
	yv1= yxf(L2v,2);
	yv2=yhat(L2v,2);  % verif period data

	nc=length(yc1);
	nv=length(yv1);

	[a1,a2,r95yc1]=acf(yc1,2);
	c(1,:)=[mean(yc1) median(yc1)  std(yc1) min(yc1) max(yc1)  a1(1) nc];

	[a1,a2,r95yc2]=acf(yc2,2);
	c(2,:)=[mean(yc2) median(yc2)  std(yc2) min(yc2) max(yc2)  a1(1) nc];

	[a1,a2,r95yv1]=acf(yv1,2);
	c(3,:)=[mean(yv1) median(yv1)  std(yv1) min(yv1) max(yv1)  a1(1) nv];

	[a1,a2,r95yv2]=acf(yv2,2);
	c(4,:)=[mean(yv2) median(yv2)  std(yv2) min(yv2) max(yv2)  a1(1) nv];

	 cc=c(1,1);  % calib-period mean of observed data
	 cv=c(3,1);  % verification-period mean of observed data



% Compute axes limits and straight lines for scatterplots
	
	aa1=min([min(yc1)  min(yc2)]);  % minimum for axes
	aa2=max([max(yc1)  max(yc2)]);  % maximum for axes
	A1 = [aa1 aa2 aa1 aa2]; % axes for scatterplot, calib period

	aa1=min([min(yv1)  min(yv2)]);  % minimum for axes
	aa2=max([max(yv1)  max(yv2)]);  % maximum for axes
	A2 = [aa1 aa2 aa1 aa2]; % axes for scatterplot, verif period
	
	[a1,a2,yl2]=lintrnd(yv2,yv1); % verification period straight line fit	

% Sign-test computations

	s1(3,:)=[0 0];  % initialize to no significance
	
	[s1(2,1),s1(1,1),a3,a4]=signtest(yv2,yv1,1);  % test on departures
	if a3=='Y', s1(3,1)=.95; end;
	if a4=='Y', s1(3,1)=.99; end;

	[s1(2,2),s1(1,2),a3,a4]=signtest(yv2,yv1,2);  % test on 1st diff
	if a3=='Y', s1(3,2)=.95; end;
	if a4=='Y', s1(3,2)=.99; end;


% Compute mean square error and related statistics

	[mae rmse re]= rederr(c(1,1),c(3,1),yc2,yc1);  % for calibration period
	s2(2,2:9)=[mae rmse re];
	[mae rmse re]= rederr(c(1,1),c(3,1),yv2,yv1);  % for verification period
	s2(1,2:9)=[mae rmse re];
	s2(:,1)=[nv nc]';  % sample size in col one

% Computations for correlation statistics
%
% F-test for significance of multiple correlation coefficient and 
% correlation coefficient from p. 93, 113 in Panofsky and Brier (1968).
% Modified so that equation for sample F and evaluation of df use
% a sample size adjusted for first-order autocorrelation (WMO 1966).

%  Adjusted size is
% n'= n * (1-r1)/ (1+r1), where n is original sample size, and r1 is the
% maximum absolute value of first order autocorrelation coefficient of
% the two series being compared. 

% Adjustment is made only if first-order autocorrelation coef of at
% least one of the two series is significant at the 95% level.  This
% restriction avoids decimating the df just because of a chance high
% autocorrelation in a small sample.

	rtemp =corrcoef(yv1,yv2);  % compute correl coef
	s3(2,1) = rtemp(1,2); 
	s3(:,2)= s3(:,1) .^2; %  R-squared and r-squared
	s3(:,7)=[nc nv]';
	
	% Adjust sample size for autocorrelation
		
	L3=[abs(c(1,6))>abs(r95yc1)  abs(c(2,6))>abs(r95yc2)];
	if any(L3);  % One of the series has significant autocorrelation
	     rc=max([abs(c(1,6))  abs(c(2,6))]);  % highest first order r for calib pd
		s3(1,3)=rc;
		ncp= nc * (1 - rc) / (1+rc)
		s3(1,8)=ncp;   % adjusted sample sizes
	else;  % no significant autocorr in either series
		s3(1,3)=0.000;
		ncp=nc;
		s3(1,8)=nc;
	end

	L4=[abs(c(3,6))>abs(r95yv1)  abs(c(4,6))>abs(r95yv2)];
	if any(L4);  % One of the series has significant autocorrelation
		rv=max([abs(c(3,6))  abs(c(4,6))]);  % highest first order r for verif pd
		nvp= nv * (1 -rv)  / (1+rv)
		s3(2,3)=rv;
		s3(2,8)=nvp;
	else;  % no signif autocorrelation
		s3(2,3)=0.000;
		nvp=nv;
		s3(2,8)=nv;
	end

	s3(1,4)= (s3(1,2) * (ncp-1-length(Ib)))  / ((1-s3(1,2)) * length(Ib));
	dfn=length(Ib);   dfd=ncp-length(Ib)-1;
	if dfd<0, error('NEGATIVE DF FOR R-SQ AFTER R1 ADJUSTMENT'); end;
	s3(1,5)=table2(f95,dfd,dfn);	

	s3(2,4)=(s3(2,2) * (nvp-2)) / (1 - s3(2,2));
	dfn=1;  dfd=nvp-2;
	if dfd<0, error('NEGATIVE DF FOR r-SQ AFTER R1 ADJUSTMENT'); end;
	s3(2,5)=table2(f95,dfd,dfn);

else
end

k2=1;

while k2~=8;
k2=menu('CHOOSE ONE','SCATTERPLOTS','HISTOGRAMS',...
  'TIME SERIES PLOTS FOR CALIB AND VERIF PERIODS',...
  'TIME SERIES PLOTS RELATED TO RE STATISTIC',...
  'SUMMARY STATISTICAL TABLES',...
  'TIME SERIES PLOTS OF RECONSTRUCTION',...
  'KEYBOARD MODE FOR PLOTS','QUIT');

if k2==1;  % scatterplots
	axis('square');
	subplot(111);
	axis(A1);
	plot(yc2,yc1,'*',[A1(1)  A1(2)], [A1(3) A1(4)],'-');
	ax=axis;  % get axes
	axis;  % unlock axes
	pp=[ax(1)+(ax(2)-ax(1))/10   ax(4)-(ax(4)-ax(3))/8];
	title('SCATTERPLOT, CALIB PERIOD');
	ylabel('OBSERVED');
	xlabel('MODEL');
	text(pp(1),pp(2),['R = ',num2str(s3(1,1))]);
	if k4=='Y', keyboard;
	else, pause, end;

	axis(A2);
	plot(yv2,yv1,'*',[A2(1)  A2(2)], [A2(3) A2(4)],'-',...
		yv2,yl2,'--');
	ax=axis;  % get axes
	axis;  % unlock axes
	pp=[ax(1)+(ax(2)-ax(1))/10   ax(4)-(ax(4)-ax(3))/8];
	title('SCATTERPLOT, VERIF PERIOD');
	ylabel('OBSERVED');
	xlabel('MODEL')
	text(pp(1),pp(2),['r = ',num2str(s3(2,1))]);
	if k4=='Y', keyboard;
	else, pause, end;
	clg
	subplot(111)
	axis
	axis('normal')
	clc; home;

elseif k2==2;  % histograms
	subplot(221)
	L=min([min(yc1) min(yv1) min(yc2)  min(yv2)]);
	M=max([max(yc1) max(yv1)  max(yc2)   max(yv2)]);
	x=linspace(L,M,10);
	[n1,xx]=hist(yc1,x);
	[n2,xx]=hist(yv1,x);
	[n3,xx]=hist(yc2,x);
	[n4,xx]=hist(yv2,x);

	ymax=max([max(n1)  max(n2)  max(n3)  max(n4)]);

	[xxx1,yy1]=bar(xx,n1);
	[xxx2,yy2]=bar(xx,n2);
	[xxx3,yy3]=bar(xx,n3);
	[xxx4,yy4]=bar(xx,n4);

	ymin=0;
	xmin=min(xxx1);
	xmax=max(xxx1);
	axis([xmin xmax ymin ymax]);

	plot(xxx1,yy1);
	title('CALIB-PERIOD OBSERVED');
	plot(xxx2,yy2);
	title('VERIF-PERIOD OBSERVED');
	plot(xxx3,yy3);
	title('CALIB-PERIOD MODEL');
	plot(xxx4,yy4);
	title('VERIF-PERIOD MODEL');
	if k4=='Y', keyboard;
	else, pause, end;
	subplot(111)
	axis
	clc; home;

elseif k2==3; % ts plots
	subplot(211)
	plot(yr1,yc2,'-',yr1,yc1,'--')
	title('CALIBRATION PERIOD:  SOLID=MODEL, DASHED=OBSERVED');
	xlabel('YEAR')

	plot(yr2,yv2,'-',yr2,yv1,'--g');
	title('VERIFICATION PERIOD: SOLID=MODEL,  DASHED=OBSERVED');
	xlabel('YEAR')
	if k4=='Y', keyboard;
	else, pause, end;
	clg
	subplot(111);
	clc; home;


elseif k2==4;  % Plots related to RE statistic
	subplot(211);
	plot(yr2,yv2,'-',yr2,yv1,'--',yr2,cc(ones(nv,1),:),'--g');
	ax=axis;  % get axes
	axis;  % unlock axes
	pp=[ax(1)+(ax(2)-ax(1))/10   ax(4)-(ax(4)-ax(3))/8];
	title(' SOLID=MODEL, DASHED=OBSERVED, HORIZ=CALIB-PD MEAN');
	xlabel('YEAR')
	text(pp(1),pp(2),['RE = ',num2str(s2(1,8))]);


	plot(yr2,yv2,'-',yr2,yv1,'--',yr2,cv(ones(nv,1),:),'--g');
	title(' SOLID=MODEL, DASHED=OBSERVED, HORIZ=VERIF-PD MEAN');
	xlabel('YEAR')
	text(pp(1),pp(2),['RE = ',num2str(s2(1,9))]);

	if k4=='Y', keyboard;
	else, pause, end;
	clg
	subplot(111);
	clc; home;


elseif k2==5;  % summary statistical tables
	clc; home; disp('BASIC STATISTICS FOR OBSERVED AND MODEL DATA');
	disp('')
	disp(...
'      MEAN    MEDIAN    ST DEV      MIN      MAX       r1        N');
	disp('');
	disp(c);
	disp(''); disp('');
	disp('     Row 1: Calibration period, observed data');
	disp('     Row 2: Calibration period, model    data');
	disp('     Row 3: Verification period, observed data');
	disp('     Row 4: Verification period, model data');
	disp('')
	disp('The above array is c, a return argument of recver');
	disp(''); disp('');
	
	disp('Any key to continue')
	pause; clc; home;

	
	disp('CORRELATION ANALYSIS:  OBSERVED WITH MODEL')
	disp('')
	disp('     R         R-sq   F-level      F-95      F-99     N-adj')
	disp('     r         r-sq                                        ')
	disp('')
	RR1=[s3(:,1)  s3(:,2)  s3(:,4)  s3(:,5)  s3(:,6) s3(:,8)];
	disp(RR1)
	disp('');
	disp('    First row above applies to calibration period.');
	disp('    Second row applies to verification period.');
	disp('');  disp('');
	disp('Any key to continue')
	pause; clc; home;

	disp('    AUTOCORRELATION ADJUSTMENT USED IN ABOVE TABLE')
	disp('')
	disp('      N     r1   N-adj')
	disp('')
	RR2=[s3(:,7)   s3(:,3)    s3(:,8)];
	disp(RR2)
	disp('') ;  disp('');
	disp('First row in above tables applies to calibration period.');
	disp('Second row applies to verification period.');
	disp('')
	disp('r1 is either 0 or the maximum of the absolute values of the');
	disp('first-order autocorrelations for the pair of series being');
	disp('compared.  Zero is used if neither first-order r is');
	disp('significant at the 95% level, as determined in acf.m');
	disp('')
	disp('Any key to continue')
	pause; clc; home;

   disp('SIGN TEST RESULTS');
	disp(' ')
	disp('DEPARTURES    1ST-DIFFERENCES');
	disp('');
	disp(s1);  
	disp(''); disp('');
	disp('Row 1 = total number of cases');
	disp('Row 2 = test statistic (p. 330, Fritts 1976)');
	disp('Row 3 = sig level: 99%, 95%, or 0 (not signif)');
	disp(''); disp('');
	disp('Any key to continue')
	pause; clc; home;

	disp('REDUCTION OF ERROR AND RELATED STATISTICS');
	disp('');
	disp(...
'    N         MAE1       MAE2      MAE3      RMSE1     RMSE2    RMSE3')
	disp('')
	disp(s2(:,1:7));   disp(''); disp('');


	disp('Row 1 = Statistics for verification period');
	disp('Row 2 = Identical statistics for calibration period');
	disp('')
	disp('MAE - MEAN ABSOLUTE ERROR')
	disp('RMSE - ROOT MEAN SQUARE ERROR')
	disp('')
	disp('   1 - for reconstruction')
	disp('   2 - if calib-pd mean of observed data were reconstruction')
	disp('   3 - if verif-pd mean of observed data were reconstruction')
	disp('')
	disp(['RE VS CALIB-PERIOD MEAN = ',num2str(s2(1,8))]);
	disp(['RE VS VERIF-PERIOD MEAN = ',num2str(s2(1,9))]);
	disp('')
disp('Any key to continue');
	pause; clc; home;

elseif k2==6;  % Time series plots of reconstruction
	recplot(yhat(:,2),s2(2,5),c(1,1),[yhat(1,1)  yhat(length(yhat),1)],k4);
	clc; clg; home;

elseif k2==7;  % keyboard mode for plots
	k4=[];
	k4=input('Keyboard-mode pauses after each plot? Y/N [N]','s');
	if isempty(k4), k4='N'; end;
else ;  % quit
	break
end
end;  % while	


