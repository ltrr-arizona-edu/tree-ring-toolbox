function [cvx,k,a,b] = cfnegx1(yrv,yrvn,xvn)
% 
% Modified on cfnegx.m to return [] as all output if cannot fit the modified negativ exponential
% 
% yrv -- years at which smoothed values are needed
% yrvn -- years at which original time series has data;  some might
%	be NaN
% xvn -- original time series;  some might be NaN
%
% cvx -- values of smooth curve at years yrv
%
% hwarn 1 if want warning window, 0 if not
%
% Fits data  to the equation
% g(t) = k + a * exp(-b*t),
%   where t is the shifted time variable t = yrvn-yrvn(1)+1
%   In other words, t is same length as yrvn after 
%   dropping NaNs
%
%
% Any NaNs in yrvn and xvn should be internal, not at the ends of the
% series

tall = yrv-yrvn(1)+1; % Shifted time variable, all years, including
 % those with NaN data

% Remove NaNs
yrvn(isnan(xvn))=[];  % delete rows with xvn as NaN
xvn(isnan(xvn))=[];

   
% Scale time series to begin at time t=1
t = yrvn-yrvn(1)+1;

% Fill data matrix, year in col 1, data in col 2
W = [t xvn];

% Estimate parameters
x0=0.01;  % Starting value of nonlinear parameter
%fopts=[0 1e-4 1e-4]; % THIS THE MOST RECENTLT USED WITH CALL TO LEASTSQ
%fopts=[0 1e-2 1e-2];;  NOT USED
% fopt.Display='off';
% fopt.TolCon=1e-4;  % Termination tolerance on the constraint violation [ positive scalar ]
% fopt.TolFun=1e-4; % Termination tolerance on the function value [ positive scalar ]
fopt=optimset('Diagnostics','off','Display','off','TolCon',1e-4,'TolFun',1e-4);

[x]=lsqnonlin('negxpk1',x0,[],[],fopt,W);
%[x] = leastsq('negxpk1',x0,fopts,[],W);
[f,g,c]=negxpk1(x,W); % f is residuals, g is growth curve, c(1 x 2)
   % is linear parameters

% Estimated parameters
k=c(1);
a=c(2);
b=x;

% OLD code from cfnegx
% if (any([k a b]<0)) & hwarn ==1,
%   parstr=['k=',num2str(k),' a=',num2str(a),' b=',num2str(b)];
%  warndlg('Neg Exp not right: ', parstr)
%  pause
%  disp('Press Any Key to Continue')
% end

% Get smoothed curve at all time values, including those that might have
% been marked off with NaNs

if (any([k a b]<0)); % if mod NE innappropriate
    cvx=[];
else;
    cvx = k + a*exp(-b*(tall));
end;

