function [f,p]=frfreq (N,pds)

% Fourier frequencies and periods in frequency range of interest
%
% D. Meko 9-29-94
%
%
%******** INPUT ARGS********
%
% N	(1 x 1) length of time series (padded length if padded)
% pds	(1 x 2) shortest period of interest, longest period of interest
%
%
%
%********** OUTPUT ARGS *********
%
% f (mf x 1) Fourier frequencies
% p (mf x 1) Corresponding periods (years)


N2 = N/2;  % There will be this many total Fourier freq, not counting 0
f1= ((1:N2)/N)';  % Fourier freqs
p1 = 1 ./ f1;  % corresp periods

L1  = p1 >= pds(1) & p1 <= pds(2); % points to subset in freq-range of interest
p=p1(L1);
f=1 ./ p;



