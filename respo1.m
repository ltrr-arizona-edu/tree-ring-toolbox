% respo1.m    response function analysis

% PCA regression equations coded from Mardia et al., 1979, p. 244-246

%*******  USER-WRITTEN FUNCTIONS AND MAT FILES NEEDED
%
% acf.m
% zscores.m
% prepmos.m
% stepr1.m
% sos.m
% sideby.m
% whit1.m
% tdist.mat --- table of student-t dist
% ftable95.mat --- table of 95% point for F-dist


%******************************************

if ~exist('tdist'), load tdist, end;  % load table for t-distrib
if ~exist('f95'), load ftable95, end; % load tbl for F-distrib


%********  INPUT DATA AND CALL PREPMOS.M TO GET LAGGED MONTHLY
%          CLIMATE ARRAY
P=input('Name of  precip array?  ','s');
T=input('Name of temperature array?  ','s');
TRI = input('Name of tree-ring array? ','s');
eval(['P= ',P,';']);
eval(['T= ',T,';']);
eval(['TRI= ',TRI,';']);
yrs1=input('Start, end years for tree array? ');
yrs4=input('Start, end years for correlation analysis? ');
endmo=input('Ending month for tree-climate year? 1=jan:  ');
nmos=input('Number of months of climate year?  ');

k4=[];
k4=input('Keyboard-mode pause for plotting acf of residuals? Y/N [N]','s');
if isempty(k4), k4='N'; end;

[A,mlabs]=prepmos(P,T,yrs4,endmo,nmos); % Organize the monthly data

[mA,nA]=size(A);
yrs2=[A(1,1)  A(mA,1)];  % Start, end year of lagged P-T array


LA=A(:,1) >= yrs4(1)  &  A(:,1) <= yrs4(2);
yrs3=(yrs1(1):yrs1(2))';  % cv of years for tree array
LTRI=yrs3 >= yrs4(1)  & yrs3 <= yrs4(2);


iw=input('Which tree-ring series to use?  ');
tree=TRI(LTRI,iw);  % tree-ring index for selected years, site

%*********  PREWHITEN TREE RINGS (OPTION) USING MODEL FIT
%      TO THE DATA ANALYSIS PERIOD COVERED BY yrs4

k3=[];
k3=input('PREWHITEN TREE-RINGS BEFORE ANALYSIS? Y/N [Y]','s');
if isempty(k3), k3='Y'; end;
if k3=='Y'
	maxar=input('MAXIMUM ORDER OF AR MODEL TO CONSIDER: ');
	iLT=find(LTRI);  % row subscripts of TRI for selected years 
	LTRI1=LTRI;   % Initialize 0/1 pointer to selected tree years with
%     AR starter-values included
	LTRI1(iLT(1)-maxar:iLT(1)-1)= ones(maxar,1); % Turn on "1s" of starter
% 			years for pointer
	ysub=TRI(LTRI1,iw);  % Tree series with starter years included.
	[eyar,arorder,vrat,arcs]=whit1(ysub,maxar,1);
	text(1,.9,['SITE ',int2str(iw),'  AR(',num2str(arorder),')']);
	text(1,.85,['res var/ orig var = ',num2str(vrat)]);

	if k4=='Y', keyboard;
	else, pause, end;
	tree=eyar(maxar+1:length(eyar));

else
end

clg, clc, home;
disp('WORKING ..DO NOT PRESS A KEY!');

%*****************  SIMPLE CORRELATION ANALYSIS, TREE VS EACH
%                   MONTHS CLIMATE

C=[tree  A(LA,2:nA)];

R1=corrcoef(C);
[mR,nR]=size(R1);
R=R1(1,2:nR);   % simple correl coefs, tree ring with various mos climate

% Compute 2-SE confid limits for correl coefs, adjusting for
% persistence.  
ry=zeros(1,nA-1);  % Preallocate
NN=zeros(1,nA-1);
cl=zeros(1,nA-1);
N=mA;  % Unadjusted sample size (number of years in A)
[r,SE2,r95]=acf(tree,2);
rx=r(1);   % r1 for tree-ring index
for i=2:nA;   % Loop for each months climate
	[r,SE2,r95]=acf(A(LA,i),2);
	ry(i-1)=max(abs([r(1) rx]));
	ryy=ry(i-1);
	NN(i-1)=((1-ryy) / (1+ryy)) * N;
	if NN>N, NN=N; end;  % Dont increase degrees of freedom for neg ar coef
	cl(i-1)=2.0 / sqrt(NN(i-1));
end



[mR,nR]=size(R);
[xx,yy]=bar(R);
ii=nA-1;


%*********   END OF CORRELATION ANALYSIS


Z=zscores(A(:,2:nA));  % standardize lagged clim array 
[mz,nz]=size(Z);
y=zscores(tree);   % Get desired tree-ring series, and standardize

%***********  MLR OF TREE VS MONTHLY CLIMATE VARIABLES

%  Stepr.m adds variables according to relative size of correlation
%   of variables not yet in model with residuals from current model. 
%   Variables entry stops at maximum adjusted R-squared.

nullfg1=0;  % change this to 1 if no significant predictors
I1=1:nz;   % Pointer to cols of climate array
[I2,stats,c,e,yhat]=stepr1(Z,y,I1); % stats(1,2,3) holds
%    R-squared, adjusted-R-squared, and F-level R-squared
%    (Draper and Smith 1981, p. 91;  Panofsky and Brier 1968, p. 113)

if(stats(1)==0), nullfg1=1; end;  % Null model
if nullfg1==0;    %  At least one predictor entered model; proceed
F1=stats(3);
dfn1 = length(I2);  dfd1=mz-length(I2)-1;  % numer and denom deg freed.
ff1=table2(f95,dfd1,dfn1);  % 95 point from f-table


beta=zeros(1,nz);  % Initialize rv for standardized regression coef
beta(I2)=c(2:length(c));
[xx3,yy3]=bar(beta);


% 95 % confidence bands around mlr coefficients (Draper and Smith 1981,
%  p. 94.  Note cautions for interpretation,  p. 95.

Zp=Z(:,I2);  % subset of Z containing only selected predictors
cvfull= (inv(Zp'*Zp)) * (std(e))^2;  % covariance matrix for reg. coefs
cvdiag = sqrt(diag(cvfull));  % estim sdevs of coeffs
t95pnt1 = table1(tdist,mz-(length(c)-1));
t95pnt2=t95pnt1(1);
c95 = t95pnt2 *  cvdiag;

else;  % nullfg1=1, relationship too weak;  no predictors enter
	clc, home, disp('WARNING: RELATIONSHIP TOO WEAK, NO PREDICTORS ENTER');
end


%**********  PCA OF MONTHLY CLIMATE;  PCA REGRESSION TO PREDICT TREE


nullfg2=0;  % this will change to 1 if no significant predictors
[RR,V,L,S,F]=eigen1(Z,1);  % pca on correl mtx of Z
W=Z*V;  %  Amplitudes of PCs of Z (cols)

% Compute variance explained of climate by each climate PC

ts=sum(diag(L));
teach=100 * diag(L) ./ ts;


a=W\y;  % Estim coefs of model to predict tree from clim eig amplitudes
ep = y - W*a;  % Compute regression residuals

% Compute pct variance of tree growth explained by each climate 
% eigenvector.

Rtemp=corrcoef([y W]);
V1a=(Rtemp(1,2:nz+1)) .^2;  % proportion variance expld
V1b= sum(V1a);  % Total proportion variance expld by model with all pc's

epsq=ep' * ep;  % sum of squares of residuals of regr of tree variable on the
			  % full set of  pc amps of climate array.

term1=epsq(ones(nz,1),:); % Terms needed for eq 8.85, p. 245 in 
aa=(mz-nz-1);             % Mardia et al. 1979
term2=aa(ones(nz,1),:);
term3= sqrt(mz * diag(L));

% Find .975 and .995 probability points of "t" distribution:  for 95%
% and 99% confidence limits around regression coefficients a.
% Only the 95% cl is used in the plots, but the 99% is there if you need
% it.

t9599=table1(tdist,aa);
t1=t9599(1);
t2=t9599(2);
t95=((t1 * term1 ./ term2) .^2) ./  term3;
t99=((t2 * term1 ./ term2) .^2) ./  term3;

[xx1,yy1]=bar(a);
xxx1=(1:nz)';   % abscissa values for plots of conf lims around a's

J = abs(a) < abs(t95);  % Ones point to elements to zero out
Jr = ~J;  % 0/1 cv, 1s point to nonzero elements of a
ev = epsq / mz;  %  Variance of residuals from PC regression using all
%		PCs as predictors.
vexp=sum(V1a(Jr));  % Total proport. variance tree growth explained by the 
%	restricted set of climate eigenvectors
Jrf = find(Jr);  % Subscripts of cols of A corresp to restricted set
%      of PC predictors

%******  F-level and significance for PC regression using only
%        coefs sig at 95% level



if sum(Jr)==0, nullfg2=1; end;  % No significant PC coefficients--null model
if nullfg2==0;    %  At least one signif pc coefficient; proceed
dfn2=sum(Jr);  dfd2=mz-sum(Jr)-1;  % deg freedom, numer and denomin
F2=vexp * (mz-sum(Jr)-1) / ((1-vexp)*(sum(Jr))) ; % F-ration for sample
ff2=table2(f95,mz-sum(Jr),sum(Jr));  % Lookup 95% point of F


%***********  TRANSFORM THE PC-REGRESSION WEIGHTS BACK INTO 
%***********  WEIGHTS ON THE INDIVIDUAL MONTHLY CLIMATE VARIABLES, AND
%***********  COMPUTE ASSOCIATED ERROR BARS (MARDIA ET AL. 1979, P. 246)

% Retain only those components whose  tree vs pc  regression coefs were
% significant at 95% level.  Form a restricted least squares estimator
% by setting all non-significant a's to zero.

arest=a;  % initialize restricted estimator
arest(J)=zeros(sum(J),1);
b=V * arest;  %   Transformed coefs -- on the original monthly climate
bv=zeros(Jr);  % Preallocate

Ld=diag(L);

for i=1:nz;  % for each monthly climate variable
	bv(i) = (ev / mz) * sum ((V(i,Jr) .^2) ./ (Ld(Jr))');
end

bv2se=2.0 * sqrt(bv);  % Two-standard errors for coefs of transformed
%		model.

dff2=[sum(Jr)-1    mz-sum(Jr)];  % degrees of freedom for f test


[xx4,yy4]=bar(b);
else;  % No significant pc-regression coefs -- null model
	clc, home;
	disp('NO SIGNIFICANT PC-REGRESSION COEFS IN MODEL');
	disp('ALL PCs TOGETHER EXPLAIN FOLLOWING PROPORTION VARNC: ');
	disp(V1b);
	pause
end;  % of code dealing with Null model



%**************   BEGIN MENU SECTION OF OPTIONAL PLOTS

k2=1;

while k2~=7;  % Stay in plot loop
clg
subplot
k2=menu('SELECT ONE','CORRELATIONS, TREE VS CLIMATE',...
    'MLR OF TREE VS MONTHLY CLIMATE',...
    'SUMMARY OF PC REGRESSION, TREE VS CLIMATE',...
	'BAR PLOTS OF INDIVIDUAL CLIMATE PCs',...
	'RESPONSE FUNCTION','KEYBOARD MODE TOGGLE','QUIT');

if k2==1;  % Correlation analysis

	plot(xx,yy,'-',1:ii,cl,'--',1:ii,-cl,'--',1:ii,zeros(1,ii),'-');
	ax=axis;   % Get coords of ends of axes
	axis;  % Unlock axes
	title('TREE VS CLIMATE CORRELATIONS, WITH 2SE BARS')
	xlabel('MONTHLY CLIMATE VARIABLE')
	ylabel('r')
	text(0.2,ax(3),['PERIOD ',int2str(yrs4(1)),' to ',int2str(yrs4(2))]);


	for i=1:nR
		if abs(R(i)) >.05
			text(i,R(i),mlabs(i,1));
		else
			text(i,0.05,mlabs(i,1));
		end
%		text(i,.050,mlabs(i,2))
%		text(i,.050,mlabs(i,3))
		text(i,.0,mlabs(i,4))
	end
	if k4=='Y', keyboard;
	else, pause, end;

elseif k2==2;  % Summary of mlr, tree vs monthly climate
	
	if nullfg1 == 0;  % At least one predictor entered model; proceed
	plot(xx3,yy3,'-',(1:nz)',zeros(nz,1),'-r',...
        I2,c95,'*',I2,-c95,'*');
	ax=axis;  % get coords of axes for plot
	axis; % unlock axes
	title('MLR, TREE RINGS VS MONTHLY CLIMATE');
	xlabel(' CLIMATE VARIABLE');
	ylabel('STDZD REGRESSION COEFFICIENT');

		for i=1:nz
		if beta(i) >.05
			text(i,beta(i),mlabs(i,1));
		else
			text(i,0.05,mlabs(i,1));
		end
	%	text(i,.050,mlabs(i,2))
	%	text(i,.050,mlabs(i,3))
		text(i,.0,mlabs(i,4))
		end
	yp= (ax(4)-ax(3))/20;  % compute y plotting pt
	
	text(0.5,ax(4)-yp,[num2str(length(I2)),' PREDICTORS EXPLAIN ',...
       num2str(stats(1)),' OF VARIANCE']);
	text(0.5,ax(4)-2*yp,['F = ',num2str(F1),'  df= ',...
	  num2str(dfn1),',',num2str(dfd1),'  F-95 = ',num2str(ff1)]);
	text(0.5,ax(3)+1*yp,['SAMPLE SIZE, N = ',int2str(mz)]);
	text(0.5,ax(3),'* -- 2-SE CONF BANDS');
	if k4=='Y', keyboard;
	else, pause, end;
	else;   % No predictors entered model;  cannot make plot
		clc, home;
		disp('NO PLOT BECAUSE NO PREDICTORS ENTERED IN MLR');
	end;   % of if dealing with null model (no predictors)

elseif k2==3;  %  PC REGRESSION, TREE VS CLIMATE
	if nullfg2 == 0;  % At least one predictor entered model; proceed
	plot(xx1,yy1,'-',xxx1,t95,'--g',xxx1,-t95,'--g',xxx1,...
            zeros(nz,1),'-r');
	title('PC REGRESSION COEFS AND 95% CONF LIMITS')
	xlabel('CLIMATE EIGENVECTOR')
	ylabel('COEFFICIENT');
	ax=axis;
	axis;  % unlock axes
	 
	yp= (ax(4)-ax(3))/20;  % compute y plotting pt
	
	text(0.5,ax(4)-yp,[num2str(sum(Jr)),' PCs EXPLAIN ',...
	  num2str(vexp),' OF VARIANCE']);
	text(0.5,ax(4)-2*yp,['F = ',num2str(F2),'   df= ',...
	  num2str(dfn2),',',num2str(dfd2),'  F-95 =  ',...
	  num2str(ff2)]);
	text(0.5,ax(3),['SAMPLE SIZE = ',int2str(mz)]);
	
	if k4=='Y', keyboard;
	else, pause, end;

	bar(V1a)
	ax=axis;
	axis;   %unlock axes
	title('TREE VARIANCE EXPLAINED BY EACH CLIMATE PC')
	xlabel('PC NUMBER');
	ylabel('PROPORTION VARIANCE');
	yp= (ax(4)-ax(3))/20;  % compute y plotting pt
	
	text(0.5,ax(4)-yp,['ALL ',num2str(nz),' PCs EXPLAIN ',...
     num2str(V1b),' OF VARIANCE']);
	text(0.5,ax(4)-2*yp,['SAMPLE SIZE = ',int2str(mz)]);
	if k4=='Y', keyboard;
	else, pause, end;

	else;   % No predictors entered model;  cannot make plot
		clc, home;
		disp('NO PLOT BECAUSE NO SIGNIF PC PREDICTORS ');
	end;   % of if dealing with null model (no predictors)

elseif k2==4;  % Plot individual climate pcs
	k5=1;
	while k5~=2;
	k5=menu('SELECT ONE','PLOT ANOTHER PC','GO BACK TO PREVIOUS MENU');
	if k5==1
		kwhich=input('WHICH PC?  ');
		VV=V(:,kwhich);
		[xx2,yy2]=bar(VV);
		plot(xx2,yy2,xxx1,zeros(nz,1),'-r');
		title(['CLIMATE PC # ',int2str(kwhich),...
		'  PCT CLIMATE VAR EXPLD= ',num2str(teach(kwhich))]);
		xlabel('MONTHLY CLIMATE VARIABLE');
		ylabel('WEIGHT')


		for i=1:nz
		if VV(i) >.05
			text(i,VV(i),mlabs(i,1));
		else
			text(i,0.05,mlabs(i,1));
		end
	%	text(i,.050,mlabs(i,2))
	%	text(i,.050,mlabs(i,3))
		text(i,.0,mlabs(i,4))
	end
		if k4=='Y', keyboard;
		else,pause, end;
	else
		k5=2
	end
	clc
	end;  % of while on k5

elseif k2==5;  % Response function plot	

if nullfg2 == 0;  % At least one predictor entered model; proceed
plot(xx4,yy4,'-',(1:nz)',bv2se,'--r',(1:nz)',-bv2se,'--r',(1:nz)',...
    zeros(nz,1),'-r');
title('RESPONSE FUNCTION, WITH 2-SE CONF. BANDS FOR COEFFICIENTS')
xlabel(' CLIMATE VARIABLE');
ylabel(' COEFFICIENT');


		for i=1:nz
		if b(i) >.05
			text(i,b(i),mlabs(i,1));
		else
			text(i,0.05,mlabs(i,1));
		end
	%	text(i,.050,mlabs(i,2))
	%	text(i,.050,mlabs(i,3))
		text(i,.0,mlabs(i,4))
		end
	if k4=='Y', keyboard;
	else, pause, end;

	else;   % No predictors entered model;  cannot make plot
		clc, home;
		disp('NO PLOT BECAUSE NO SIGNIF PC PREDICTORS ');
	end;   % of if dealing with null model (no predictors)

elseif k2==6;  % Keyboard-mode toggle to make easy to get hardcopy plots
	k4=[];  % Next 2 lines allow execution to be interrupted after each
%    screen plot.  Facilitates making  .met files.
	k4=input('Keyboard-mode pauses after each plot? Y/N [N]','s');
	if isempty(k4), k4='N'; end;


elseif k2==7;  % Quit the plot loop
	break
end
end;  % of while on k2

