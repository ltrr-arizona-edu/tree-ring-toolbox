function [b,r2,yhat] = regress1(y,X)
%
% Modified from MATLAB function regress.m by D Meko on 8-28-95
%
%REGRESS Performs multiple linear regression using least squares.
%	(X is an nxp matrix, y is the nx1 vector of observations.) 
%	References:
%	   [1] Samprit Chatterjee and Ali S. Hadi, "Influential Observations,
%	   High Leverage Points, and Outliers in Linear Regression",
%	   Statistical Science 1986 Vol. 1 No. 3 pp. 379-416. 
%	   [2] N. Draper and H. Smith, "Applied Regression Analysis, Second
%	   Edition", Wiley, 1981.

%	B.A. Jones 3-04-93
%	Copyright (c) 1993 by The MathWorks, Inc.
%	$Revision: 1.4 $  $Date: 1993/10/04 12:26:29 $

if  nargin < 2,              
    error('REGRESS requires at least two input arguments.');      
end 


% Check that matrix (X) and left hand side (y) have compatible dimensions
[n,p] = size(X);
[n1,collhs] = size(y);
if n~=n1, 
    error('The number of rows in Y must equal the number of rows in X.'); 
end 

if collhs ~= 1, 
    error('Y must be a vector, not a matrix'); 
end

% Find the least squares solution.
[Q R]=qr(X);
b = R\Q'*y;

yhat=X*b;

% Calculate R-squared.
RSS=norm(yhat-mean(y))^2;
TSS=norm(y-mean(y))^2;
r2=RSS/TSS;
end


