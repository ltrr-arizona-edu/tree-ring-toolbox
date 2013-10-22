function [L1,I1]=getnear1(NY,W,npref)
%
% Given a matrix NY of number of modeling overlap years for
% each of m1 chronologies with n1 climate stations, 
%
[m1,n1]=size(NY);

% compute threshold number of required years.  
nthresh=min(max(NY'));

% Compute row index to chrons that have satisfy npref years
L1=NY>=npref;
I1 =  find((any(L1')')
