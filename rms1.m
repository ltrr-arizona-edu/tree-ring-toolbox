function s=rms1(y,yh,p)

% Residual mean square (Weisberg 1985, 12)

% Say you've already run a regression, and got the predicted values.
% The residual mean square can be used to put confidence intervals
%   around the predicted values, as in Weisberg, p. 229.


%*************************   INPUT ARGS
%
% y (N x 1)  observed values of predictand for a calibration period
% yh (Nh x 1) predicted values for the calibration period.
% p (1 x 1)  number of parameters in model (includes constant term)
%



%*********************  OUTPUT ARGS *********
%
% s  (1 x 1) residual mean square, defined as the sos of residuals
%     divided by the df.

if (nargin ~= 3), error ('Need 3 input args'), end;

N=length(y);
N1=length(yh);

if N ~= N1
	error ('Length of y must equal length of yh')
end



% Compute residuals.

e= y - yh;  % residual equals obs minus predicted

s= sum(e .* e) / (N-p);
