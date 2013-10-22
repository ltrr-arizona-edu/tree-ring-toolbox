function [TR,TI,Rtilde,varchk]=pdpb(x,p,k,padlen)
%
% Fast Fourier transform using Peter Bloomfield method. Fortran
% program on p 113-114 of: 
%
% Peter Bloomfield, 1976. "Fourier Analysis of
% Time Seris, An Introduction", John Wiley & Sons, Inc., 258 pp.
%   p 113-114 -- Fortran pgm for periodogram.
%
%******************  IN ARGS *************************
%
% x (N x 1)r  time series, length N years (or whatever time units)
% p (1 x 1)r  total proportion of time series to be tapered
% k (1 x 1)i  detrending option: 0=subtract mean, 1=subtract least-
%		squared fit straight line
% padlen (1 x 1)i desired padded length of x; should be a power of 2
%		and be at least as large as N, but no larger than 32X the 
%		next largest power of 2
%
%*************** OUT ARGS ******************************
%
% TR (? x 1)r Real part of Fourier transform
% TI (? x 1)r Imaginary part Fourier transform
% Rtilde (? x 1)r  Blomfield's squared amplitude of Fourier transform,
%		which is proportional to the periodogram.  Note that first
%		np2/2+1 elements are returned 
% varchk (1 x 2)r  variance computed from time series and from
%		periodogram
%
%
%************* USER-WRITTEN FUNCTIONS CALLED *******************
%
% taper -- taper ends of time series
%
%
%**************  METHOD ***************************************
%
% -Detrend by subtracting mean or least-squares fit straight line
% -Optionally taper the detrended time series
% -Compute the padded length, np2, as the next power of 2 equal to or
%	longer than the original series length, N
% -Pad the tapered series to length np2 by appending zeros
% -Compute the FFT of the padded, tapered seriesd
% -Compute the squared amplitude from the FFT. From the eq on bottom of 
%  p. 113, would multiply squared amplitude by n/(8*pi) to get 
%  true periodogram, I(omega). 
% -Check that properly scaled sum of squared amplitudes equals variance
%  of the  time series (after tapering, but before padding)
%
%********************* NOTES ****************************************
%
% Relation between squared amplitudes, Rtilde, and periodograms.
%
% Bloomfield notes that Rtilde multiplied by N/(8*pi) gives the
%		Bloomfield periodogram
%
% The frequency points for Rtilde start at 0 and have increment of
%	1/np2 in "per year" units or 2*pi/np2 in radians
%
% The Splus spec.pgram periodogram is obtained by multiplying
%		Rtilde by N/4.  It follows that the Splus periodogram is
%		2*pi times the Bloomfield periodogram
%
% Calling function should check that padlen valid:  at least as large
%		as original series length, a power of 2, and no larger than
%		32X the next power of 2 above the series length

% Size and check arguments
[mx,nx]=size(x);
if mx<20 | nx ~=1,
	error('x must be col vector of minimum length 20')
end
if k<0 | k>1,
	error('Detrending parameter k must be 0 or 1')
end
if p<0 | p>=1; % taper proportion is capital P in Bloomfield
	error('Taper proportion must satisfy 0< p <=1')
end	

N = mx; % Time series length; to match Bloomfield variable name


% Detrend by subtracting mean or least-squared fit straight line
x=dtrend(x,k);

% Taper the time series, if desired
if p~=0,
	x = taper(x,p);
end


np2=padlen;  % just to fit in with Bloomfield notation





% Compute FFT
xx= zeros(np2,1); % Next two statements amount to padding with zeros
xx(1:N)=x;
J = fft(xx,np2); % Z is complex matrix, fourier transform of
		% xx

X=real(J);
Y=imag(J);

npgm = np2/2+1; % number of periodogram ordinates
con = (2.0 * np2/N) ^2;

% Compute Bloomfield's 'squared amplitude'
Rtilde = (J .* conj(J)) * con; % Equivalent to (X .^2 + Y .^2)* con


% Adjust transforms and Rtilde because Bloomfield's notation has a factor 
% of 1/n in summations for FFT, while matlab function fft.m does not
TR=X/np2;
TI=Y/np2;
Rtilde = Rtilde / (np2*np2); % Adjusment due to matlab fft function not
		% scaling by 1/np2 in summations for FFT, while Bloomfield
		% notation as factor 1/np2 in summations

% Because of symmetry in Rtilde, only need first np2/2+1 values
Rtilde=Rtilde(1:npgm); % first "half" of squared amplitudes
TR=TR(1:npgm);
TI=TI(1:npgm);

% Variance check.  See variance equalities, p. 50 in Bloomfield.
% Compute variance as sum of squares of departures from mean divided 
% by N. Compare with variance computed from squared amplitudes. If 
% tapering was used, the variance is that ofthe tapered series.
var1 = var(x)*(N-1)/N; % Adjustment because built-in variance function
	% for matlab uses N-1 divisor
var2 =(1/2)*(N/np2)*(sum(Rtilde(2:npgm-1))+Rtilde(1)/2+Rtilde(npgm)/2);
varchk=[var1 var2];


