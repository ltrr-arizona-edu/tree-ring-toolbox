function [mae,rmse,re]=rederr(ybarc,ybarv,yhat,y)
% rederr:  verification statistics for a climate reconstruction
% [mae,rmse,re]=rederr(ybarc,ybarv,yhat,y);
% Last revised 4-22-92
%
%*** INPUT ARGS
%
% ybarc (1 x 1r  calibration period mean of observed predictand
% ybarv (1 x 1)r  verification-period mean of observed predictand
% yhat (m1 x 1)r  verification-period time series of reconstructed predictand
% y (m2 x 1)r verification-period time series of observed predictand
%
%*** OUTPUT ARGS
%
% mae (1 x 3) mean absolute error assuming the reconstruction is:
%      mae(1):  the reconstructed time series itself (no assuming here)
%      mae(2):  ybarc: the calibration-period mean of the observed data
%      mae(3):  ybarv: the verification-period mean of the observed data
% rmse (1 x 3) root-mean-square error, ordered as in mae
% re   (1 x 2) reduction-of-error statistic
%     re(1) --  based on mean square errors for the reconstruction and
%			compared with mse if reconstruction were a constant equal
%			to ybarc
%     re(2) --  based on mean square errors for the reconstruction and
%			compared with mse if reconstruction were a constant equal
%			to ybarv
%
%
%*** REFERENCES
% 
% The mae and rmse error are defined in many standard time series and regression texts, 
% for example:
% Weisberg, S., 1985, Applied Linear Regression, 2nd ed., John Wiley, New York.
%
% The re statistic in a tree-ring context is discussed in:
% Gordon, G., 1982, Verification of dendroclimatic reconstructions, In:  
% M.K. Hughes et al. (eds.), Climate from tree rings.  Cambridge University Press,
% Cambridge, UK, 58-61.
%
%*** UW FUNCTIONS CALLED
%*** TOOLBOXES NEEDED
%
%*** NOTES
%

%*****     Size arrays and do a few consistency checks

if nargin~=4  &  nargout~=3
	error('Incorrect number of input or output arguments')
end

[m1,n1]=size(yhat);
[m2,n2]=size(y);
L1 = [n1~=1  n2~=1 m1~=m2];
if any(L1)
	error('yhat, y must be same-length CVs')
end

n=length(y);  % sample size


% Compute errors of reconstruction, mean abs errors,mean square err
% and root-mean-square error

e1=yhat-y;   % reconstruction error
e2=ybarc-y;  % errors if ybarc were reconstruction
e3=ybarv-y;  % errors if ybarv were reconstruction

mae=[mean(abs(e1))  mean(abs(e2))  mean(abs(e3))];
mserr=[e1'*e1   e2'*e2  e3'*e3] / n;  % mean square errors
rmse = sqrt(mserr);

% Compute re statistics

re(1) = (mserr(2) - mserr(1))/ mserr(2);
re(2) = (mserr(3) - mserr(1))/ mserr(3);
