function p= cohprob1(gsq)
% p= cohprob1(gsq)
%
% 95% and 99% probability points for distribution of squared coherency
%
% D Meko 7-22-95
%
%*************** IN ARGS ***************************************
%
% gsq (1 x 1)r  sum-of-squares of filter used to smooth periodograms
%		and cross-periodogram
%
%
%*************** OUT ARG ****************************************
%
% p (1 x 3)r  90%, 95% and 99% probability points of distribution of 
%		squared coherency
%
%
%***************** METHOD *************************************
%
% Bloomfield, Peter, 1976. Fourier Analysis of Time Series:  An 
% Introduction.  John Wiley & Sons, 258 pp -- page 227

a = gsq/(1-gsq); % exponent for equation preceding eq 13 in 
	% Bloomfiel, p. 227.

p= 1 - (1-[.90 .95 .99]) .^a
