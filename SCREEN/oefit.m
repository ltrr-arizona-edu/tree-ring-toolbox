function [s,r2,ob,of,f1,f2,way,k1,smore]=oefit(Z,yr,nn,nranked,kopt)
% [s,r2,ob,of,f1,f2,way,k1,smore]=oefit(Z,yr,nn,nranked,kopt);
%
% Fit output-error (OE) model to fulle period of data. Unlike oesplit.m,
% oefit does not use split-sample modeling to choose best order model.
% Strategy is to pick model that has lowest final prediction error, with
% constraint that all parameters at least 2 std devs from zero. If no parameters
% sig different from zero, take OE(1,0) by default
%
% D Meko, 12-5-97
%
%************************* IN ARGS **************************
%
% Z (mZ x 2)r time series of output (col 1) and input (col 2)
% yr (myr x 1)i year vector for Z
% nn (? x 3)i orders nb,nf,nk for each candidate model
% nranked (1 x 1)i number of 'best' models to consider
% kopt(1 x 1)i  1=regular mode,  2=diagnostics mode
%
%************************* OUT ARGS *************************
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
%		chosen.  Automatically []!  The following applied to oesplit.m
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
% k1 (1 x 1)i  indicator whether all parameters of final model 
%     significant (more than 2 std devs from zero); 1==yes, 0==no
% smore (1 x 1)r  incremental gain in explained variance above that
%     of the OE(1,0) model.  NaN if the chosen model is OE(1,0)
%
%********************** NOTES **************************
%
% Call "A" models those calibrated on the first half and validated
% on the second half; and "B" models those calibrated on the
% second half and validated on the first half.
%
% If some model has all parmeters significant,
% the order of preference is:
%	1-ranking according to minimum FPE
%	2-if a tie by "1", highest validation variance-explained fraction
%
% Split sample results in equal samples if total number of 
% years is even.  If total number is odd, first half is 1 larger
% than the second

a=NaN;

% Check, size and allocate
[m1,n1]=size(Z);
if m1~=length(yr) | n1~=2,
	error('Z mus be 2-col matrix with same row size as yr')
end
nyrs=length(yr);


[m1,n1]=size(nn);
if m1<2 | n1~=3
	error('nn must have minimum of 2 rows and must have 3 cols')
end
nmods=m1;  % number of candidate OE models

% Reduce nranked, number of ranked models to be considered, if
% that number as an input argument is too large for nn
if nmods<nranked
	nranked=m1;
end



%***************************  FIT MODELS

k1=0;

% Subtract the  column means
mn = mean(Z);
z = Z - repmat(mn,nyrs,1);;
vary = var(z(:,1)); % variance of the output



%-----------  OE(1,0) MODEL

nnn=[1 0 0]; % specify model order

% Estimate the model parameters
th = oe(z,nnn);
thmod1=th;

% Calculate flag for whether **any** significant parameters,
% where significant is defined as two or more sdevs from 
% zero
% Get variances of estimate parameters; convert to standard
% deviations; test for greater than two standard deviations
% from zero
nparams=th(1,5)+th(1,8); % number of B params + number of F params
cp = th(3,1:nparams); % parameter estimates, B-operator, then F
vcp = th(4:(3+nparams),1:nparams); % variances of parameters
c2 =2*sqrt(diag(vcp)); % two standard deviations of parameters
k1mod1= any(abs(cp)>abs(c2'));
k2mod1 = all(abs(cp)>abs(c2')); % all params signif

% Variance computation
ysim = idsim(z(:,2),th); % noise-free simulation
yh=ysim;  % this is also y-hat
vare =  var(z(:,1)-ysim); % variance of the errors
smod1 = 1 - (vare/vary); % Explained variance stat
fpemod1=th(2,1); % final prediction error

% Correlation of predicted with actual
rr2 = corrcoef([z(:,1) yh]);
r2mod1 = rr2(1,2);


%-----------  OE (2,0) MODEL

nnn=[2 0 0]; % specify model order

% Estimate the model parameters
th = oe(z,nnn);
thmod2=th;

% Calculate flag for whether **any** significant parameters,
% where significant is defined as two or more sdevs from  zero
% Get variances of estimate parameters; convert to standard
% deviations; test for greater than two standard deviations
% from zero
nparams=th(1,5)+th(1,8); % number of B params + number of F params
cp = th(3,1:nparams); % parameter estimates, B-operator, then F
vcp = th(4:(3+nparams),1:nparams); % variances of parameters
c2 =2*sqrt(diag(vcp)); % two standard deviations of parameters
k1mod2= any(abs(cp)>abs(c2'));
k2mod2 = all(abs(cp)>abs(c2')); % all params signif

% Variance computation
ysim = idsim(z(:,2),th); % noise-free simulation
yh=ysim;  % this is also y-hat
vare =  var(z(:,1)-ysim); % variance of the errors
smod2 = 1 - (vare/vary); % Explained variance stat
fpemod2=th(2,1); % final prediction error

% Correlation of predicted with actual
rr2 = corrcoef([z(:,1) yh]);
r2mod2 = rr2(1,2);



%-------------  OE (1,1)

nnn=[1 1 0]; % specify model order

% Estimate the model parameters
th = oe(z,nnn);
thmod3=th;

% Calculate flag for whether **any** significant parameters,
% where significant is defined as two or more sdevs from  zero
% Get variances of estimate parameters; convert to standard
% deviations; test for greater than two standard deviations
% from zero
nparams=th(1,5)+th(1,8); % number of B params + number of F params
cp = th(3,1:nparams); % parameter estimates, B-operator, then F
vcp = th(4:(3+nparams),1:nparams); % variances of parameters
c2 =2*sqrt(diag(vcp)); % two standard deviations of parameters
k1mod3= any(abs(cp)>abs(c2')); % any signif params
k2mod3 = all(abs(cp)>abs(c2')); % all params signif

% Variance computation
ysim = idsim(z(:,2),th); % noise-free simulation
yh=ysim;  % this is also y-hat
vare =  var(z(:,1)-ysim); % variance of the errors
smod3 = 1 - (vare/vary); % Explained variance stat
fpemod3=th(2,1); % final prediction error

% Correlation of predicted with actual
rr2 = corrcoef([z(:,1) yh]);
r2mod3 = rr2(1,2);

%*********************  PICK THE "BEST" MODEL

fpe = [fpemod1 fpemod2 fpemod3];
[xmin,imin]=min(fpe);

if imin==1; % OE(1,0);
   bestmod=1;
   k1=k1mod1;
	r2=r2mod1;
	s=smod1;
	fpe=fpemod1;
	ob=1;
	of=0;
   thmod=thmod1;
elseif imin==2; %OE(2,0)
   bestmod=2;
   k1=k1mod2;
	r2=r2mod2;
	s=smod2;
	fpe=fpemod2;
	ob=2;
	of=0;
	thmod=thmod2;
elseif imin==3; % OE(1,1)
   bestmod=3;
   k1=k1mod3;
	r2=r2mod3;
	s=smod3;
	fpe=fpemod3;
	ob=1;
	of=1;
   thmod=thmod3;
end


way=[];


% Store info on supplemental variance explained above OE(1,0) model
if bestmod==1;
   smore = NaN;
else
   smore = s - smod1;
end




% Compute number of significant (99% level) autocorrelations
% of residuals and of crosscorrelations between input and
% residuals.  Consider only lags 1-10
[e,r,f1,f2]=resid2(z,thmod,11);



%***********  optional detailed-look info 

if kopt==2;
	file1=uiputfile('*.mat','Save file for detailed look');
	set1=[' z mn yr  s  smore ob of f1 f2 way k1 '];
	eval(['save ',file1,set1])
end

