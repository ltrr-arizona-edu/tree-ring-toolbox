function H=hatmtx(X)
% hatmtx: hat matrix from multiple linear regression
% CALL: H=hatmtx(X);
%
% Meko 3-29-98
%
%*************  IN
%
% X (mX x nX)r  matrix of predictors, including col of ones as first col
%    if a constant term in the model
%
%**************** OUT 
%
% H (mH x nH)r  hat matrix
%
%*****************  NOTES
%
% Source:  Weisberg, 1985, p. 109


H = X * (inv(X' * X)) * X';
