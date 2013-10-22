function [sss,epps]=wigley1(r,NN)
% wigley1:  subsample signal strength and expressed population signal
% [sss,epps]=wigley1(r,NN);
% Last revised 9-2-99
%
% In tree-ring standardization via treei, computes SSS and EPS
%
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
% epps (NN x 1)r  expressed population signal for chronologies of
%		sample size 1,NN
%
%
% REFERENCES
% Wigley et al. (1984)
%
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%***  NOTES  **********************************
%
% In tree-ring application, sample size is number of trees, not
% number of cores.
%___________________________________________________________


% Compute SSS
n = (1:NN)';
N = NN(ones(NN,1),:);
top = n .* (1 + (N-1)*r);
bottom = N .* (1 + (n-1)*r);
sss = top ./ bottom;


% Compute EPS
top = n*r;
bottom =  1+ (n-1)*r;
epps = top ./ bottom;


% End of file
