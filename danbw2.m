function [bw,df]=danbw2(w,nobs,npad);
% danbw2: bandwidth and degrees of Daniell filter for smoothing periodogram into spectrum
% CALL:  [bw,df]=danbw2(w,nobs,npad);
%
% Meko  3-27-98
%
%************************ IN 
%
% w (nw x 1)r  weights of the Daniell filter
% nobs (1 x 1)i  number of observations in the time series spectrum computed for
% npad (1 x 1)i  padded length of the time series (as in pdgm5.m) spectrum computed for
%
%*********************** OUT 
%
% bw (1 x 1)r  bandwidth, in "per year" frequency units 
% df (1 x 1)r  df  degrees of freedom for resulting spectral estimates
%
%*********************************  NOTES
%
% If a series has length nobs, the Fourier frequencies are spaced at 1/nobs.
% If a series has been padded with zeros to length npad, the F frequencies are at 1/npad
% Approach is to find the rectangular filter with same sum of squares as the 
% Daniell filter, and use the bandwidth and df of that rectangular filter as the
% bandwidth and df of the Daniell filter.  The computed df and bw are used in context
% of spectral analyis for smoothing the periodogram computed, say, in a calling 
% function.  In computing the bw and df, must adjust for possible padding of the
% time series 

%--------  compute sum of squares of the Daniell filter
sos1 = sum(w .* w);

n5=1/sos1; % Length of rect filter with same sum of squares as the Daniell filter
% call this the equivalent rectangular filter

% Bandwidth of the equavalent rectangular filter is simply the number of weights
% of the rectangular filter times the "Fourier" frequency of the time series.
% If the series has been padded with zeros, this Fourier frequency must be computed 
% using the padded length rather than the original length
bw = n5/npad;  % Bandwidth of the rectangular filter with the same sum of squares as
% the input Daniell filter
% The degrees of freedom of a rectangular filter applied to the periodogram of an
  
% unpadded time series of length nobs is  two times the number of weights divided by
% the series length.  The df must be adjusted down if the series was padded
% before computing the periodogram
df = (n5*2)* (nobs/npad);
