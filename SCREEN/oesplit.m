function [s,r2,ob,of,f1,f2,way,k1]=oesplit(Z,yr,nn,nranked,kopt)
% [s,r2,ob,of,f1,f2,way]=oesplit(Z,yr,nn,nranked,kopt);
%
% Fit output-error (OE) model by split-sample method, selecting 
% best model as the one with minimum error variance in validation.
%
% D Meko, 4-27-96
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

% Initialize some storage matrices; The columns of SA and SB will
% be as follows:
% 1 row-index to model structure matrix nn
% 2 final prediction error (a calibration statistic)
% 3 validation variance-explained fraction
% 4 "all parameters significant" flag
% 5 total number of B-order and F-order parameters
% 6 total number of B-order parameters
SA=a(ones(nmods,1),ones(6,1)); % general stats, "A" models
SB=a(ones(nmods,1),ones(6,1)); % general stats, "B" models
FA=a(ones(nmods,1),:); % Final prediction error, "A" models
FB=a(ones(nmods,1),:); % Final prediction error, "B" models


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

% ********** CALIBRATE ON FIRST HALF, VALIDATE ON SECOND HALF
%

% Pull out the calibration and validation subsets
ze=Z(LCA,:);
zv=Z(LVA,:);

% Subtract the calibration-period mean
mn = mean(ze);
MNe = mn(ones(nyrsA,1),:);
MNv = mn(ones(nyrsB,1),:);
ze = ze - MNe;
zv = zv - MNv;

% Compute variance of the validation-part of the output series.
% Will need this later for variance-explained statistic
vary = var(zv(:,1)); 


% Loop over candidate OE models
for n = 1:nmods;
	kmod = nn(n,:);
	th = oe(ze,kmod);
	[yh,fit]=compare(zv,th,1);
	FA(n)=th(2,1);
	SA(n,2)=FA(n);
	SA(n,3)=1-((fit*fit)/vary);
	% Store total number of parameters, and total number ofB-order
	% parameters
	nparams=th(1,5)+th(1,8);
	SA(n,5)=nparams;
	SA(n,6)=th(1,5);
	% Get variances of estimate parameters; convert to standard
	% deviations; test for greater than two standard deviations
	% from zero
	cp = th(3,1:nparams); % parameter estimates, B-operator, then F
	H = th(4:(3+nparams),1:nparams); % variances of parameters
	c2 =2*sqrt(diag(H)); % two standard deviations of parameters
	SA(n,4)= all(abs(cp)>abs(c2'));
end

% Rank models according to final prediction error
% and pull out nranked top-ranked models
[FA,iv]=sort(FA);
iv = iv(1:nranked);
SA = SA(iv,:); % re-order the statistics matrix
SA(:,1)=iv; % row index to model-structure matrix nn


% ********** CALIBRATE ON SECOND HALF, VALIDATE ON FIRST HALF
%
% Pull out the calibration and validation subsets
ze=Z(LCB,:);
zv=Z(LVB,:);

% Subtract the calibration-period mean
mn = mean(ze);
MNe = mn(ones(nyrsB,1),:);
MNv = mn(ones(nyrsA,1),:);
ze = ze - MNe;
zv = zv - MNv;

% Compute variance of the validation-part of the output series.
% Will need this later for variance-explained statistic
vary = var(zv(:,1)); 

% Loop over candidate OE models
for n = 1:nmods;
	kmod = nn(n,:);
	th = oe(ze,kmod);
	[yh,fit]=compare(zv,th,1);
	FB(n)=th(2,1);
	SB(n,2)=FB(n);
	SB(n,3)=1-((fit*fit)/vary);
	% Store total number of parameters, and total number ofB-order
	% parameters
	nparams=th(1,5)+th(1,8);
	SB(n,5)=nparams;
	SB(n,6)=th(1,5);
	% Get variances of estimate parameters; convert to standard
	% deviations; test for greater than two standard deviations
	% from zero
	cp = th(3,1:nparams); % parameter estimates, B-operator, then F
	H = th(4:(3+nparams),1:nparams); % variances of parameters
	c2 =2*sqrt(diag(H)); % two standard deviations of parameters
	SB(n,4)= all(abs(cp)>abs(c2'));
end

% Rank models according to final prediction error
% and pull out nranked top-ranked models
[FB,iv]=sort(FB);
iv = iv(1:nranked);
SB = SB(iv,:); % re-order the statistics matrix
SB(:,1)=iv; % row index to model-structure matrix nn


%***************** PICK OVERALL BEST MODEL ****************

% Store information on top-ranked models
nA=SA(1,5); % total number params
nAA=SA(1,6); % number of B-operator params
nB=SB(1,5);
nBB=SB(1,6);
iA=SA(:,1); % row index to model-structure matrix nn
iB=SB(:,1); % ditto

% Information on whether any model has all params significant
LA=SA(:,4)==0;
LB=SB(:,4)==0;
L=all(LA) & all(LB); % no model has all params significant

% If no model has all parameters significant, pick the simplest
% top-ranked model
if L; % if no model has all parameters significant
	if nA<nB; % fewer total params
		nnn=nn(iA(1),:); % selected model
	elseif nB<nA,
		nnn=nn(iB(1),:);
	else
		if nAA<=nBB,
			nnn=nn(iA(1),:);
		else
			nnn=nn(iB(1),:);
		end
	end
	way=6;
else; % At least one "A" or "B" model has all params significant
	% Change LA and LB so that a "1" means all parameters signif
	LA=~LA;
	LB=~LB;
	if any(LA) & any(LB); % At least one "B" model and one "A"
		jA = min(find(LA));
		jB=min(find(LB));
		s1=SA(jA,3); % validation variance-expld fraction
		s2=SB(jB,3);
		nnnA=nn(SA(jA,1),:);
		nnnB=nn(SB(jB,1),:);
		if s1>=s2
			nnn=nnnA;
			if jA==1,
				way=2;
			else
				way=3;
			end
		else
			nnn=nnnB;
			if jB==1,
				way=2;
			else
				way=3;
			end
		end
		if jA==jB & jA==1 & all(nnnA==nnnB);
			way=1;
		end
	elseif any(LA) & ~any(LB);  % An "A" has all signif, no "B" does
		jA = min(find(LA));
		nnn=nn(SA(jA,1),:);
		if jA==1,
			way=4;
		else
			way=5;
		end
	elseif any(LB) & ~any(LA);  % A "B" has all signif, but no "A" does
		jB = min(find(LB));
		nnn=nn(SB(jB,1),:);
		if jB==1,
			way=4;
		else
			way=5;
		end
	end
	

end



%************ FINISHED WITH SPLIT-SAMPLE MODELING.  NOW FIT
% SELECTED BEST MODEL TO WHOLE DATA

k1=0;

% Subtract the means
mn = mean(Z);
MNz = mn(ones(nyrs,1),:);
z = Z - MNz;

% Estimate the model parameters
th = oe(z,nnn);

% Calculate flag for whether **any** significant parameters,
% where significant is defined as two or more sdevs from 
% zero
% Get variances of estimate parameters; convert to standard
% deviations; test for greater than two standard deviations
% from zero
nparams=th(1,5)+th(1,8);
cp = th(3,1:nparams); % parameter estimates, B-operator, then F
vcp = th(4:(3+nparams),1:nparams); % variances of parameters
c2 =2*sqrt(diag(vcp)); % two standard deviations of parameters
k1= any(abs(cp)>abs(c2'));

% Variance computation
%vare = th(1,1); % estimate variance of residuals -- questionable result
ysim = idsim(z(:,2),th); % noise-free simulation
vare =  var(z(:,1)-ysim); % variance of the errors
vary = var(z(:,1)); % variance of the output
s = 1 - (vare/vary); 

% Compute no-noise prediction
[yh,fit] = compare(z,th,1);

% Correlation of predicted with actual
rr2 = corrcoef([z(:,1) yh]);
r2 = rr2(1,2);


% Store the orders of the "B" and "F" parameters of the 
% model
ob = th(1,5); 
of = th(1,8);

% Compute number of significant (99% level) autocorrelations
% of residuals and of crosscorrelations between input and
% residuals.  Consider only lags 1-10
save hocus z th 
[e,r,f1,f2]=resid2(z,th,11);
