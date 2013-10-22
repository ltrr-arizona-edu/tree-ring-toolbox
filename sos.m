function [c,stats]=sos(y,X);
% sos: sum of squares terms for multiple linear regression
% CALL: [c,stats]=sos(y,X);
% 
% Meko 4-18-97
%
%*********************  IN **********************
%
% y (n x 1)r dependent variable
% X (n x m1) independent variables
%
%
%********************  OUT *********************
% 
% c ((p+1) x 1)r regression coefs, constant first
% stats (1 x 3)r R-squared, adjusted R-squared, F-ratio
%
%***************** USER FUNCTIONS CALLED *************
%
% stepr.m stepwise regression
%
%**************** NOTES **************************
% 
% Reference: Draper and Smith 1981, p. 89-96


[n,m1] = size(X);
p=m1+1;  % total number of parameters includes a constant

n=length(y);  % number of observations
tss=y' * y;  % total sum of squares
ssb0=length(y) * (mean(y))^2;   % Sum of squares due constant term
ctss=tss-ssb0;

c1=ones(length(y),1);
X1=  [c1 X];
c = X1\y;  % regression coefs, constant first

ssr= c' * X1' * y;  % regression ss
ssrb0= ssr - ssb0;  % regression ss given constant term b0
ssres = tss - ssr;  % residual ss


fact = ssres ./ ctss;
rsq=  1 - fact;  %  r-squared statistic
rsqadj = 1 - fact * (n-1)/(n-p);  % adjusted r-squared

Frat= (ssrb0/(p-1)) / (ssres/(n-p));  % F-ratio for equation

stats=[rsq rsqadj Frat];

