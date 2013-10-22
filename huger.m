function [G,a,b,g]=huger(w);

% Hugershoff function

n=length(w);  % length of input time series
t=(1:n)';
X= [ones(n,1)  log(t)  t];  % predictor matrix
y=log(w);   % predictand vector

c= X\y;  % estimate coefs
yhat=X*c;  % predicted natural log series

G=exp(yhat);  % Computed growth curve

% Transform coefficients
a=exp(c(1));
b=c(2);
g=c(3);
