function [f,g,c]=negexpk(x,W,nm,xlab,ylab)

% Nonlinear least squares fit of modified negative exponential curve
%
% D. Meko 12-6-93
%
% Fits data in second column of W to the equation
% g(t) = k + a * exp(-b*t)
% Also plots the original data, the fitted curve, and annotates plot
% 	with series label, axes labels, parameter estimates and error norm.
% 
% Two linear parameters and one nonlinear parameter.
%  c(1) = k
%  c(2) = a
%  x(1) = b
%
%
%******************   INPUT ARGS   ******************
%
% x (1 x 1) initial value of nonlinear parameter
% W (mW x 2) year (or other nominal x-axis variable) in col 1
%	Data to be fit in col 2.
% nm (1 x ?) str name of series -- used to label plot
% xlab (1 x ?) str xlabel for plot -- e.g., year
% ylab (1 x ?) str ylabel for plot -- e.g., Ring Width (mm)
%
%
%*******************  OUTPUT ARGS ******************
%
% f (mW x 1) residuals = observed minus fitted
% g (mW x 1) fitted curve
%
%
%*********  USE ***********************************
%
%  Generally used with leastsq.m.  Form of call is :
%    x=leastsq('negexpk',x0,[],[],W,nm,xlab,ylab)
%    x0 is the initial value of x for starting iterations.
%  Also could be called as stand-alone after using leastsq.m to get
%    nonlinear parameter estimate.  In this way, second call could use
%    last nonlinear estimate as the input parameter value x0.


[mW,nW]=size(W);
yr=W(:,1);
w=W(:,2);


t=(1:mW)';
A=ones(mW,2);  % initialize predictor matrix

A(:,2)=exp(-x(1)*t);  % exponential term

c=A\w;  % c(1) will be constant term;  c(2) coefficient of exp term
g=A*c;   % fitted growth curve
f=w - g;  % residuals = observed - fitted

plot(yr,w,yr,g);
xlabel(xlab);
ylabel(ylab);
xt2 = yr(1)+(yr(mW)-yr(1))/5;;
xt1 = yr(1)+(yr(mW)-yr(1))/6;
yt =6* max(w)/8;
text(xt2,1.3*yt,['g(t)= ',num2str(c(2)),'  exp(-',num2str(x(1)),' t) + ',num2str(c(1))])
text(xt2,1.2*yt,['err norm = ' sprintf('%g',norm(f))])
text(xt1,1.4*yt,nm);
