function [hstar,H,sep,extrap] =sepred2(X,L3,sereg)
% sepred2:  standard error of prediction and extrapolation indicator for regression
% CALL:  [hstar,H,sep,extrap] =sepred2(X,L3,sereg);
%
% Meko 3-30-98
%
%******************  IN 
%
% X (ms x ns)r predictor matrix, including calib pd and earlier data (ones column as col 1)
% L3 (mL3 x1)L  logical pointer to calibration period data in X
% sereg(1 x 1)r  standard error of regression
%
%****************** OUT
%
% hstar distance measure (p. 237)
% H    hat matrix
% sep  standard error of prediction (p. 229)
% extrap ....logical indicator of extrapolation (1) vs interpolation (0)
%
%*****************  NOTES
%
% Source: Weisberg (1985)
%

% Size
[mX,nX]=size(X);
yr1 = (1:mX)'; % relative "year"vector
yrs1 = [1 mX];

% Check that no missing data in calibration period
Xc = X(L3,:);
if     any(any(isnan(Xc')));
   error('NaN in calibration data');
end

% Hat matrix
H = hatmtx(X(L3,:));

% Need to merge long-term and calib data for call to mce1
yr2 = yr1(L3);
yrs2= [min(yr2) max(yr2)];


%Compute hstar for each year (Weisberg 1985, 229)
[S,hmax]= mce1(X,yrs1,yrs2);
hstar=S(:,2);

% Extrapoplations
extrap=S(:,3);

% Compute sepred as in Weisberg, 1985, eqn 9.2, p. 229
sep= sereg * sqrt(1+hstar);


