function [C,P]=crspec1(x,y,padlen,bw,spans,p,k)
%
% [C,P]=crspec1(x,y,padlen,bw)
%
% Coherency and phase for two time series
%
% D Meko 7-21-95
%
%
%******************  IN ARGS **************************************
%
% x (mx x 1)r first time series
% y (my x 1)r second time series
% padlen (1 x 1)i padded lengths of time series
% bw (1 x 1)r desired bandwidth of estimate
% spans (1 x ns)i   spans for Daniell filter
% p (1 x 1)r  total proportion of time series to be tapered
% k (1 x 1)i  detrending option: 0=subtract mean, 1=subtract least-
%		squared fit straight line
%
%************************  OUT ARGS ********************************
%
% C (mC x 1)r coherency of cross-spectrum (Bloomfield, p. 214, second eqn)
% P (mP x 1)r phase of cross-spectrum
%
%
%*************** USER-WRITTEN FUNCTIONS CALLED  *****************************
%
% pdpb -- compute FFTs and periodograms (Peter Bloomfield method)
% specpb -- spectral estimate by Peter Bloomfield method
% moddan -- smooth series with modified daniell filter
% polar -- convert real and imag transform parts to ampl and phase
%
%*****************  NOTES ***************************************
%
% 	C is the coherency, not the squared coherency!  What you will want
%  to plot for showing results is the squared coherency, C .* C


[mx,nx]=size(x);
[my,ny]=size(y);
if nx~=1 | ny~=1 | mx~=my,
	error('x and y must be col vectors of same length')
end
nobs=mx;  % number of observations in the time series before padding

% Check that padlen is valid.  Padded length must be a power of two,
% must be at least as long as the original series length, and must
% be no more than 32 times the series length.  The "32" has no special
% meaning. Who in their right mind would want to pad a 1000 year series
% to length greater than 32*1024 anyway?
if padlen<nobs, error('padlen shorter than original length'), end
ncheck = nextpow2(nobs);
next2 = 2 ^ncheck;
if padlen>32*next2, 
	error('padlen more than 32X next power of two greater than nobs')
end
n2 = next2*[1 2 4 8 16 32];
if ~any(n2==padlen)
	error('padlen not a power of 2')
end
np2=padlen;	

% Compute the Fourier transforms and periodograms of x and y
[TRx,TIx,Rtildex,varchkx]=pdpb(x,p,k,padlen); 
[TRy,TIy,Rtildey,varchky]=pdpb(y,p,k,padlen);
N = length(Rtildex);  % number of periodogram ordinates

% Compute estimated autospectra of x and y
[gx,start,step ] = specpb(Rtildex,nobs,padlen,spans,1);
[gy,start,step ] = specpb(Rtildey,nobs,padlen,spans,1);

% Scaling factor
con = (np2*np2)/ (2 * pi * nobs);


% Compute cross-periodogram
CR = (TRx .* TRy) + (TIx .* TIy); % real part
CI = (TIx .* TRy) - (TRx .* TIy); % imaginary part
TR = CR*con;
TI = CI*con;

% Smooth the cross-periodogram to get cross-spectrum
TRs = moddan(TR,1,spans);
TIs = moddan(TI,-1,spans);

% Get amplitude and phase of cross spectrun
[mag,P]=polar(TRs,TIs);

% Compute coherency
C = mag ./ (sqrt(gx .* gy));




