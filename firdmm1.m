function [b,b1,w1,A,pv]=firdmm1(pp25,nw,nfreqs)

% Hamming-windowed lowpass filter design-- filtfilt version
% D. Meko 11-3-92
% See p. 2-58 in signal processing toolbox

% Design a low-pass filter and compute and plot
% its frequency response.  User specifies the 50% freq response
% for a filter b1.  Time series will be filtered with b1 forwards and
% then backwards.  This double filtering has effect of squaring the freq
% response, so that the effective filter is b=conv(b1,b1).
% So the final filter b will have a lower frequency response at pp25 than
% 0.5.  If the filter length is long enough,  
% the period of the 25% response will approach pp25. 

%********* OTHER USER-WRITTEN FUNCTIONS NEEDED 
%
% respoint.m  - interpolates period of specified freq response
%
%
%********* INPUT ARGUMENTS
%
% pp25 (1 x 1) desired period in yrs of 50% freq response of filter
%  		b1, which will be convoluted with itself to get filter b
% nw (1 x 1) specified number of weights in filter, before convolution
%   The final filter will have 2* nw - 1 weights.
% nfreqs (1 x 1) number of frequency points between 0 and
%   0.5 per year that freq resp is to be returned for

%***********  OUTPUT ARGUMENTS
%
% b (1 x 2nw)  the filter weights convoluted withself
% b1 (1 x nw) the filter weights
% w1 (nfreqs x 1) frequencies (per year) for A
% A  (nfreqs x 1) frequency response
% pv (1x4) period (yr) of 90%, 50%, 25% and 10% freq response
%


n=nw-1;  
Wn=2.0 / pp25;  % Cutoff frequency on scale such that 0 is 
	% 0 per year, and Wn=1 is 0.5 per year (Nyquist).  At
	% this frequency, response of filter b1 is 0.5, and
	% response of b=conv(b1,b1) is 0.25 (if filter is broad 
	% enough)

b1=fir1(n,Wn);

% Compute and plot magnitude frequency response

b=conv(b1,b1);
[h,w]=freqz(b,1,nfreqs);
A=abs(h) ; % Amplitude
%ph=phase(h');  % Phase
mn=ones(nfreqs,1) * 0.5;  % horiz line at fre resp of 1/2


w1=w / (2*pi);  % want freq in per year

p90=respoint(A,w1,0.9);
p50=respoint(A,w1,0.5);
p25=respoint(A,w1,0.25);
p10=respoints(A,w1,0.1);
pv=[p90 p50 p25 p10];

plot(w1,A,w1,mn); % plot magnitude response
title('MAGNITUDE OF FREQUENCY RESPONSE')
xlabel('Frequency (per year)')
ylabel('Response')
text(0.2,0.9,['Spec. Period of 25% Response = ',num2str(pp25),' yr']);
text(0.2,0.85,['Spec. Number of Weights = ',int2str(nw)]);
text(0.2,0.8,['90% Response: ',num2str(p90),' yr']);
text(0.2,0.75,['50% Response: ',num2str(p50),' yr']);
text(0.2,0.70,['25% Response: ',num2str(p25),' yr']);
text(0.2,0.65,['10% Response: ',num2str(p10),' yr']);

pause
