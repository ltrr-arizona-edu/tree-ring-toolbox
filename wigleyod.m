%function [sss eps]=wigley1(r,NN)
% renamed wigleyod.m 2-2-96  by dmm because sohel had a more recent
% versio nelsewher
% Subsample signal strength and expressed population signal
% Source: Wigley et al. (1984)
%
%************* IN ARGS *********************************
%
% r (1 x 1) -- mean between-tree correlation coefficient for pairs
%		of core indices
% NN (1 x 1) -- full sample size on which correlation based
%
%
%*****************  OUT ARGS *******************************
%
% sss (NN x 1)r  subsample signal strength for chronologies of
%		sample size 1,NN
% eps (NN x 1)r  expressed population signal for chronologies of
%		sample size 1,NN
%
%
%********************  NOTES  **********************************
%
% In tree-ring application, sample size is number of trees, not
% number of cores.


% Compute SSS
n = (1:NN)';
N = NN(ones(NN,1),:);
top = n .* (1 + (N-1)*r);
bottom = N .* (1 + (n-1)*r);
sss = top ./ bottom;


% Compute EPS
top = n*r;
bottom =  1+ (n-1)*r;
eps = top ./ bottom;



