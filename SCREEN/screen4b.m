% script file screen4b.m
%
% Diagnostic checking of Output-Error (OE) model fits to a
% pair of input and output series as followup to running of screen3.m
% Screen3.m called oefit.m, which used FPE on full period of data
% to pick best moel 
%
% Note that screen4 (or 4a) used split-sample calib/validation to
% pick best model, while 4b uses fit to full data
% 
% If you set O(1)==2 in running screen3.m and ran the model 
% on a particular chronology, you got a .mat file with these 
% variables, which must be in the workspace to run screen4b.m:
%
% z (mz x 2)r y(t) in col 1, u(t) in col 2;  means already removed
% mn (1 x 2)r  the means for z
% yr (mz x 1)i year vector for z
% s (1 x 1)r fraction-of-variance-explained statistic
% smore(1 x1)r additional explained variance (above OE(1,0) model)
% r2 (1 x 1)r correlation between predicted and actual output
% ob (1 x 1)i order of B parameter of selected model
% of (1 x 1)i order of F parameter of selected model
% f1 (1 x 1)i number of significant (99% level) autocorrelation
%		of residual in lags 1-10
% f2 (1 x 1)i number of significant (99% level) crosscorrealation
%		between input and residuals, in lags -10 to +10
% way (1 x 1)i : indicator (1,2, or 3) for how best model was
%		chosen: as for now, only NaN is valid
% k1 (1 x 1)i  indicator whether all parameters of final
%		significant ( at least 2 stdevs from zero)
%		k1==1 yes, 
%       ==0 no

close all
nyrs=length(yr); % total number of years, calibration+validation
a=NaN;


% Input and output series, as departures from means mn
u1 = z(:,2);
y1 = z(:,1);


%*******************************************************************
% TSP OF Z-SCORES OF INPUT AND OUTPUT, FULL PERIOD
figure(1)
subplot(1,1,1)
r4=corrcoef([u1 y1]);
r4=r4(1,2);
txt4a=['r= ',num2str(r4)];
plot(yr,zscores(y1),yr,zscores(u1),'--',...
	[yr(1) yr(length(yr))],[0 0 ])
legend('Output','Input')
ylabel('Z-Scores')
title(['FULL-PERIOD INPUT AND OUTPUT; ',txt4a])



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

% Recall that mn holds full period means of output and input

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

plot(yr,y2+mn(1),yr,yh+mn(1),'--',...
[yr(1) yr(nyrs)],[mn(1) mn(1)],'w-')
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