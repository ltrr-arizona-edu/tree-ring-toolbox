% pdgm.m
% an M-file to estimate spectrum by smoothed periodogram
load F:\matlab\A.dat
nsize=282 % number of years of time series
% N is the padded length

% Load and plot the time series.
year=A(1:nsize,[1]);
x=A(1:nsize,[4]);
plot(year,x)
title('raw time series')

% subtract mean from x and compute its  discrete Fourier transform as
% defined in Ljung.  Pad to length 512 with zero's first.
% Bloomfield's DFT is 1/sqrt(n) times Ljung's DFT.
x=dtrend(x);
pause
plot(year,x)
grid

y=fft(x,512);
title('time series, mean subtracted')
grid
t=0:1:512;
f=t/512;
pause
plot(f(1:256),y(1:256))
title('discrete fourier transform, Ljung')


% Compute the periodogram as defined by Ljung by squaring the DFT
% Bloomfield's squared dft is 1/N times Ljung's DFT.
% Bloomfield's periodogram is 1/2*pi times Ljung's periodogram.
Pyy = y.*conj(y)/512;
pause
plot(f(1:256),Pyy(1:256));
title('unsmoothed Ljung periodogram')


pause
% compute the unsmoothed Bloomfield periodogram as 1/2pi * Pyy
% PyyB= (1.0/(2.0*pi)) .*Pyy;
% plot(f(1:256),PyyB(1:256));
% title('unsmoothed Bloomfield periodogram')

pause
% compute variance of time series and compare to
% 1/nsize times the sum of the Ljung periodogram ordinates
size(x)
variance = (std(x)*std(x))
chkvar = (1/(nsize-1)) * sum(Pyy(1:512))

