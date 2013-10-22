% script file oesplit2
%
% Diagnostic checking of Output-Error (OE) model fits to a
% pair of input and output series as followup to running of screen.m
% Screen.m called oesplit1.m, which used split-sample calibration-
% validation to get statistics on various models and pick a
% best model.  Screen.m operated on many series at a time. If 
% you set kopt==2 in running screen.m and restricted the model to
% run on a particular chronology, you got a .mat file with these 
% variables, which must be in the workspace to run oesplit2.m:
%
% nn (? x 3)i orders nf,nb,nk of candidate models, as from struc3.m
% Z (mZ x 2)r y(t) in col 1, u(t) in col 2
% yr (mZ x 1)i year vector for Z
% yrA (? x 1)i year vector for first "half" of Z
% yrB (? x 1)i year vector for second "half" of Z
% SA (? x 6)r, SB(? x 6)r  -- statistics from split-sample
%		calibration/validation.  "A" models were estimated
%		on first half of data and validated on second half. "B"
%		models were estimated on second half and validated on first
%		half. The rows of SA and SB were sorted and ranked
%		previously in increasing order of FPE for calibration.  Thus
%		row 1 holds the model with lowest FPE, etc.
%		The cols of SA and SB are :
%
% 	1 row-index to model structure matrix nn
% 	2 final prediction error (a calibration statistic)
% 	3 validation variance-explained fraction
% 	4 "all parameters significant" flag
% 	5 total number of B-order and F-order parameters
% 	6 total number of B-order parameters
%
% FA, FB
%
% s (1 x 1)r fraction-of-variance-explained statistic
% r2 (1 x 1)r correlation between predicted and actual output
% ob (1 x 1)i order of B parameter of selected model
% of (1 x 1)i order of F parameter of selected model
% f1 (1 x 1)i number of significant (99% level) autocorrelation
%		of residual in lags 1-10
% f2 (1 x 1)i number of significant (99% level) crosscorrealation
%		between input and residuals, in lags -10 to +10
% way (1 x 1)i indicator (1,2, or 3) for how best model was
%		chosen:
%		1-top-ranked "A" or "B" models both have all params signif;
%			and both top-ranked models are same order
%		2-top-ranked "A" or "B" model with all params significant;
%			secondary screening by validation variance-expld fraction
%		3-a lower-ranking model with all params significant; secondary
%			screening with validation variance-expld fraction
%		4-a top-ranked "A"  or "B" model; all parms signif
%			no need for screening because if "A", no significant
%			models in "B", and vice versa
%		5-lower-ranking "A" or "B" model; all params significant;
%			no need for screening ...
%		6-No model with all params significant; pick simplest top-ranked
%			model;  priority for screening: (1) lowest total number
%			of parameters, (2) lowest number of "B-operator"
%			parameters
% k1 (1 x 1)i  indicator whether any of parameters of final
%		(not split sample) model has any significant paramters
%		k1==1 yes, at least one param is signif (2 sdevs from zero)
%		k1==0 none of params is signif
%

close all
nyrs=length(yr); % total number of years, calibration+validation
a=NaN;

[m1,n1]=size(nn);
nmods=m1;  % number of candidate OE models


% Get year vectors for split samples.  The "A" split is calibrating
% on the first half of the data and verifying on the second. The
% "B" split is calibrating on the second half, verifying on the
% first.
nyrsA = ceil(nyrs/2);
nyrsB = nyrs -nyrsA;
yrA=yr(1:nyrsA);
yrB=yr((nyrsA+1):nyrs);
LCA = yr>=yrA(1) & yr <= yrA(length(yrA)); % pointer calibration-A
		% rows of Z 
LVA = yr>=yrB(1) & yr <= yrB(length(yrB));% pointer to validation-A
	%		rows in Z
LCB = LVA; % pointer to validation-B rows in Z
LVB = LCA; % pointer to calibration-B rows in Z

% Pull out the calibration and validation subsets for "A" models
ze1=Z(LCA,:); % estimation data for est on first, valid on second
zv1=Z(LVA,:); % validation data

% Subtract the calibration-period mean; get resulting "adjusted"
% series that have the first-half mean subtracted
mn1 = mean(ze1);
MNe1 = mn1(ones(nyrsA,1),:);
MNv1 = mn1(ones(nyrsB,1),:);
zze1 = ze1 - MNe1;
zzv1 = zv1 - MNv1;

% Pull out the calibration and validation subsets for "B" models
ze2=Z(LCB,:);
zv2=Z(LVB,:);

% Subtract the calibration-period mean; get resulting "adjusted"
% series that have the second-half mean subtracted
mn2 = mean(ze2);
MNe2 = mn2(ones(nyrsB,1),:);
MNv2 = mn2(ones(nyrsA,1),:);
zze2 = ze2 - MNe2;
zzv2 = zv2 - MNv2;

% Compute variance of the validation-part of the output series.
% Will need this later for variance-explained statistic
vary1 = var(zv1(:,1)); 
vary2 = var(zv2(:,1));

% Input and output series, original units
u1 = Z(:,2);
y1 = Z(:,1);
mn3=mean(Z); % full period means of output, input


%*******************************************************************
% TSP OF Z-SCORES OF INPUT AND OUTPUT, FULL PERIOD
figure
subplot(1,1,1)
r4=corrcoef([u1 y1]);
r4=r4(1,2);
txt4a=['r= ',num2str(r4)];
plot(yr,zscores(y1),yr,zscores(u1),...
	[yr(1) yr(length(yr))],[0 0 ])
legend('Output','Input')
ylabel('Z-Scores')
title(['FULL-PERIOD INPUT AND OUTPUT; ',txt4a])


%****************************************************************
% TSP OF OUTPUT (TOP) AND INPUT (BOTTOM) IN ORIGINAL UNITS, WITH
% HORIZONTAL LINES SHOWING MEANS FOR FIRST AND SECOND HALVES OF DATA
%

figure
txt1a=['A,B means = ',num2str(MNe1(1,1)),', ',num2str(MNe2(1,1))];
txt1b=['A,B means = ',num2str(MNe1(1,2)),', ',num2str(MNe2(1,2))];

subplot(2,1,1)
plot(yr,y1,yrA,MNe1(:,1),yrB,MNe2(:,1),'c--')
ylabel('Tree-Ring Index')
title(['RAW DATA AND SPLIT-SAMPLE MEANS :',txt1a]);

subplot(2,1,2)
plot(yr,u1,yrA,MNe1(:,2),yrB,MNe2(:,2),'c--')
ylabel('Precip (% of Normal)')
title(['INPUT :',txt1b]);


%*********************************************************************
% GRAPH SHOWING VARIATION OF FPE (CALIBRATION) AND VARIANCE EXPLAINED
% PERCENTAGE (VALIDATION) WITH MODEL FOR SPLIT-SAMPLE CALIBRATION
figure
txt2a=['OE(1,0)';'OE(1,1)';'OE(2,0)'];
subplot(2,1,1)
plot(SA(:,1),SA(:,2),'+',SB(:,1),SB(:,2),'x',...
   SA(1,1),SA(1,2),'yo',SB(1,1),SB(1,2),'mo');
legend('y+','First-Half Calibration','mx',...
		'Second Half Calibration','o','Lowest-FPE Model')
title(' FINAL PREDICTION ERROR IN SPLIT-SAMPLE CALIBRATION')
ylabel('FPE')
set(gca,'Xlim',[0.8 3.2],'Xtick',[1 2 3],'XtickLabels',txt2a)


subplot(2,1,2)
plot(SA(:,1),SA(:,3),'+',SB(:,1),SB(:,3),'x')
legend('y+','First-Half Calibration','mx',...
		'Second Half Calibration')
title(' 1 - VAR(ERRORS)/VAR(OUTPUT)')
ylabel('VAR-EXPLD STATISTIC')
set(gca,'Xlim',[0.8 3.2],'Xtick',[1 2 3],'XtickLabels',txt2a)

%*********************************************************************
% GRAPH SHOWING WHICH SPLIT-SAMPLE CALIBRATION MODELS HAD
% ALL PARAMETERS SIGNIFICANT
figure
txt3a=['NO ';'YES'];
subplot(2,1,1)
plot(SA(:,1),SA(:,4),'*')
title(' ALL  MODEL PARAMETERS SIGNIFICANT? -- A MODELS')
set(gca,'Xlim',[0.8 3.2],'Xtick',[1 2 3],...
	'XTickLabels',txt2a,'Ylim',[-0.25 1.25],...
	'Ytick',[0 1],'YtickLabels',txt3a)

subplot(2,1,2)
plot(SB(:,1),SB(:,4),'*')
title(' ALL  MODEL PARAMETERS SIGNIFICANT? -- B MODELS')
set(gca,'Xlim',[0.8 3.2],'Xtick',[1 2 3],...
	'XTickLabels',txt2a,'Ylim',[-0.25 1.25],...
	'Ytick',[0 1],'YtickLabels',txt3a)



%*******************************************************************
% SOME BOOKKEEPING FOR FULL-PERIOD MODELS
%
% Recall that ob and of give orders of B and F operators of 
% "best" full-period model
%
% Subtract means from input and output
u2=detrend(u1,0);
y2=detrend(y1,0);
z2=[y2 u2];

%****************************************************************
% FIT FULL-PERIOD MODEL

% Recall that mn3 holds full period means of output and input

% Fit the full-period model
th=oe(z2,[ob of 0]);


%*****************************************************************
% IMPULSE AND STEP RESPONSE FUNCTIONS

% FIRST THE EMPIRICAL IMPULSE RESPONSE
 
figure
M=20; % want for 20 lags -- lags 0 thru 19
AN=10; % prewhiten input with AR(10) model and filter output
	% by same in computing IRF
M1=M-1; % maximum lag for plotting
[IR,R,CL] = CRA(z2,M,AN,0); 
%	The "0" means no automatic plot in the above call to cra
%  Returned are the irf (IR), covariances (R), and 99% 
%  confidence level for IR
%
% Note:  you could get crosscorrelation function of y2,u2 if
% replace z2 with zscores of y2,u2;  and set AN==0 in the above
% call

stem((0:M1)',IR)
hold on
plot((0:M1)',0+CL(ones(M,1),:),'y--',...
	(0:M1)',0-CL(ones(M,1),:),'y--');

% NOW THE MODEL IMPULSE RESPONSE

ufake=zeros(M,1);
ufake(1)=1;
yfake=idsim(ufake,th)
plot((0:M1)',yfake,'m*')
title('EMPIRICAL (o) AND MODEL (*) IMPULSE RESPONSE FUNCTIONS')
text(5,CL-.005,'99% CONFIDENCE LIMIT')
set(gca,'XLim',[-0.3 20])

hold off

%******************************************************
% PLOTS OF ACF OF RESIDUALS AND CROSSCORRELATION OF INPUT WITH
% RESIDUALS, ALONG WITH 99% CONFIDENCE LIMITS
figure
resin(z2,th,25)


%****************************************************************
% AMPLITUDE AND PHASE OF FREQUENCY RESPONSE,
% WITH 2-STDEV ERROR BARS; FIRST FOR MODEL

% Get transfer function and noise model
[G,NSP] = th2ff(th);

% Pull frequency, ampli, phase, standard deviation of amp and phase
% from G
[W,AMP,PHAS,SD_AMP,SD_PHAS] = getff(G);

% Convert units of W from radians to cycles per year
W = W / (2*pi);

% Compute empirical transfer function
ge = spa(z2);
% col 1 will be frequency in radians, col 2 is transf function,
% col 3 is std dev of ampl , col 4 is phase, col 5 is std of phase
ge(1,:)=[];  % cut off row 1, the data-type indicators
We=ge(:,1)/(2*pi);


% Plot amplitude and phase of transfer function, model and 
% empirical.  Include 2-sdev bars around empirical

figure
subplot(2,1,1)

% Use this commented out version if want 2-sd devs around model
% rather than empirical tf
%plot(W,AMP,W,AMP+2*SD_AMP,'m--',W,AMP-2*SD_AMP,'m--',...
%	We,ge(:,2))

plot(W,AMP,'m',...
 We,ge(:,2),'y-',...
 We,ge(:,2)+2.0*ge(:,3),'y--',...
 We,ge(:,2)-2.0*ge(:,3),'y--')


title('FREQUENCY RESPONSE AND 2-SDEV BARS')
ylabel('AMPLITUDE')




subplot(2,1,2)
%plot(W,PHAS,W,PHAS+2*SD_PHAS,'m--',W,PHAS-2*SD_PHAS,'m--')
plot(W,PHAS,'M',...
 We,ge(:,4),'y',...
 We,ge(:,4)+2.0*ge(:,5),'y--',...
 We,ge(:,4)-2.0*ge(:,5),'y--')

ylabel('PHASE')
xlabel('FREQUENCY (1/YR)')




%*****************************************************************
% ZEROS AND POLES OF FINAL MODEL, WITH 2-SDEV ERRORS
% This plot only if ob+of>1
% In other words, zpplot will bomb for OE(1,0) model
figure
subplot(1,1,1)
if (ob+of)>1,
	zepo=th2zp(th);
	zpplot(zepo,3)
	title('Zeros (o) and Poles (x) of Model')
else
	title('No zero-pole plot possible for this model')
end


%*****************************************************************
%  NOISELESS SIMULATION OF OUTPUT FROM INPUT

figure
yh=idsim(u2,th);

% Compute variance-explained statistic
eey = y2-yh;  % model residuals
ves = 1.0 - (var(eey)/var(y2));

% correlation
temp1=corrcoef([yh y2]);
temp1=temp1(1,2);

plot(yr,y2+mn3(1),yr,yh+mn3(1),[yr(1) yr(nyrs)],[mn3(1) mn3(1)],'w-')
legend('Observed','Model Predicted')
xlabel('Year')
ylabel('Tree-Ring Index')
title('ABILITY OF MODEL TO REPRODUCE TREE-RING INDEX')

% Make a plotting point 5% down from top, 30% in from left
xtemp=get(gca,'XLim')
ytemp=get(gca,'Ylim')
x1=xtemp(1) + 0.30*(abs(diff(xtemp)));
y1=ytemp(2) - 0.05*(abs(diff(ytemp)));
txt6=['EV = ',num2str(ves)];
text(x1,y1,txt6)
x1=xtemp(1) + 0.30*(abs(diff(xtemp)));
y1=ytemp(2) - 0.10*(abs(diff(ytemp)));
txt7=['r = ',num2str(temp1)];
text(x1,y1,txt7)
