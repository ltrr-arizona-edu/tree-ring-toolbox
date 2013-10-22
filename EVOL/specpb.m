function [y,start,step ] = specpb(Rtilde,nobs,padlen,spans,sym)
%
% Spectral estimation by smoothing a periodogram
%
% D. Meko 7-19-95
%
%
%********************* IN ARGS **********************************
%
% Rtilde (mx x 1)r Bloomfield's squared amplitude (p. 48), which
%	is (R*R)/4 in 2nd eqn on p. 48, and is proportional to the
%	true periodogram (p. 48, p. 113)
% nobs (1 x 1)i  Number of original (before padding) observations
%	in the time series that Rtilde was computed for
% padlen (1 x 1)i padded length of the time series for which
%		Rtilde was computed
% spans (1 x ns)i  spans of the modified Daniell filters to be
%	used to smooth the periodogram
% sym (1 x 1)i type of symmetry to be used in filtering
%	1=even, -1=odd, 0=tack on zeros (see moddan.m)
%
%*********************** OUT ARGS ***************************
%
% y (my x 1)r estimated spectrum, same length as Rtilde
% start (1 x 1)r frequency at y(1)
% step (1 x 1)r increment of freqency in y
%
%************************ USER-WRITTEN FUNCTIONS CALLED **********
%
% moddan.m  -- modified Daniell filtering
%
% 
%************** METHOD ************************************
%
%- set a constant mult factor that to convert Rtilde to true periodogrm
%- convert Rtilde to true periodogram
%- smooth periodogram with specified modified Daniel filter
%- compute the frequency points for the spectral estimates
%
%
%*********** NOTES ***********************************************
%
% padlen is passed as an imput argument as a formality.  The padded 
% length is related to the length of Rtilde by 
% 	length(Rtilde) == padlen/2 + 1
% so padlen could have been computed in this function
%
% step is in radian units.  divide by 2*pi to get "per-year" units.
% 
% Relationship of estimated spectra to the estimated spectra from 
% SPLUS spec.pgram.  Multiply y by 2*pi, and spectra are identical.

[mx,nx]=size(Rtilde);

npgm = mx; % number of periodogram ordinates

con = nobs/(8.0*pi); % mult factor
x = Rtilde*con; % true periodogram
y=moddan(x,sym,spans); % spectral estimate

start = 0.0; % frequency for y(1);
step = 2.0*pi/padlen; % frequency increment

