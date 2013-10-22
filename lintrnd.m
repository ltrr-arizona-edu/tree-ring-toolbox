function [b,r,yhat]=lintrnd(x,y)

% Given two time series, return the linear regression equation,
% correlation coefficient, and predicted values from linear 
% regression


m1=length(x);
A=[ones(m1,1), x];
B=y;
b= A\B;  % estimate regr coefs

yhat = A * b;

r=corrcoef(x,y);
r=r(1,2);
