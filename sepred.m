% sepred.m

% Standard error of prediction from regression, plus related quantities.
% This file is identical to wrr15.m.  I have not yet made the function
% general

% Make table showing 2-se bar around min, median, max values in 
% long-term log-10 reconstruction.  And show corresp +- around
% estimate in untransformed (acre-ft) units

% tnew is the key table to be output
%   row 1 minimum,  row 2 mdian, row3 maximum
%      col 1: year
% 	    col 2: recon log-10 discharge
%      col 3: se of prediction for this log-10 discharge
%		 col 4: lower CI of antilog of col 2
%		 col 5:  antilog of col 2
%		 col 6: upper CI of antilog of col 2


% Estimated standard error of prediction for reconstruction 
% Weisberg (1985, 229, eqn 9.2)
% In discharge units; also +- bars in untranformed units (acre-ft)

% Also get:
%    s .... residual mean square
%	  sest ... standard error of the estimate (p. 12)
%	  se ... standard error of prediction (p. 229)
%	  hstar ... the distance measure (p. 237)
%	  hmax ... the max value of hstar for the calib period (p. 237)


% Assumes yhat.mat loaded

% User written functions needed:  mce1.m,   rms1.m,  pltext.m  tdist.mat
load tdist; % t-dist table

Nc=70;  % Number of years in calibration period
p=9+1; % number of parameters in model
df=Nc-p;  % 70 obs in calib pd, p  parameters (9 predictors + a constant)

%Compute hstar for each year (Weisberg 1985, 229)
[S,hmax]= mce1(X,[1663 1985],[1916 1985]);
hstar=S(:,2);


% Compute residual mean square for the log-10 reconstruction
s=rms1(log10(ycal(:,1)),log10(ycal(:,2)),p); % 10 refers to 9 predictor
  % vbls plus constant term.. Ten in df for residuals.
sest = sqrt(s); % standard error of the estimate

% Compute sepred as in Weisberg, 1985, eqn 9.2, p. 229
se= sest * sqrt(1+hstar);

yrr=(1663:1985)';

% Plot log-10 recon with 95% conf interval
% First compute multiplier from t(60) distribution
tmult=table1(tdist,Nc-p); % give both .05 and .01 levels
tmult=tmult(1);  % since want only .05 prob

plot(yrr,ylog(:,2),yrr,ylog(:,2)+tmult*se, yrr,ylog(:,2)-tmult*se);
title('Reconstructed variable and 2-se confidence bar')
ylabel('Log-10 acre-ft')
pause


% Compute tables of values of log-10 and untransformed reconstruction
% with 2-se confidence band for (1) minimum recon value, 
% (2) median, and (3) maximum

Imin=find(min(ylog(:,2))==ylog(:,2));
Imed=find(median(ylog(:,2))==ylog(:,2));
Imax=find(max(ylog(:,2))==ylog(:,2));

t15log=[ylog(Imin,1) ylog(Imin,2)  se(Imin)  ylog(Imin,2)-tmult*se(Imin) ...
      ylog(Imin,2)+tmult*se(Imin);
      ylog(Imed,1) ylog(Imed,2)  se(Imed)  ylog(Imed,2)-tmult*se(Imed) ...
      ylog(Imed,2)+tmult*se(Imed);
      ylog(Imax,1) ylog(Imax,2)  se(Imax)  ylog(Imax,2)-tmult*se(Imax) ...
      ylog(Imax,2)+tmult*se(Imax)];

tnew=[t15log(:,[1 2 3])  exp(log(10)*t15log(:,[4 2 5]))]
