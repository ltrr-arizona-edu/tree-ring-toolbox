function [w1,A,pv]=freqres1(b1,nfreqs)
% freqres1:  Hamming-windowed lowpass filter design
%
% Compute and plot frequency response of a fir filter.
% User specifies the weights in the filter and number of frequencies for the plot
% 
% D. Meko 12-15-92
%
%********* INPUT ARGUMENTS
%
% b1 (1 x nw) the filter, assumed symmetric, with nw weights
% nfreqs (1 x 1) number of frequency points between 0 and
%   0.5 per year that freq resp is to be returned for
%
%***********  OUTPUT ARGUMENTS
%
% w1 (nfreqs x 1) frequencies (per year) for A
% A  (nfreqs x 1) frequency response
% pv (1x3) period (yr) of 90%, 50% and 10% freq response
%
%********* USER-WRITTEN FUNCTIONS NEEDED 
%
% respoint.m  - interpolates period for specified frequency response
%

nw = length(b1); % number of weights in filter

% Compute and plot magnitude frequency response
b=b1;
[h,w]=freqz(b,1,nfreqs);
A=abs(h) ; % Amplitude
%ph=phase(h');  % Phase
mn=ones(nfreqs,1) * 0.5;  % horiz line at fre resp of 1/2


w1=w / (2*pi);  % want freq in per year

% Get periods for specified amplitude of frequency response
p90=respoint(A,w1,0.9);
p50=respoint(A,w1,0.5);
p10=respoint(A,w1,0.1);
pv=[p90 p50 910];

hp1=plot(w1,A,w1,mn); % plot magnitude response
cord=get(gca,'ColorOrder');
set(hp1(1),'Color',cord(1,:));
set(hp1(2),'Color',cord(1,:),'LineStyle',':');
title('AMPLITUDE OF FREQUENCY RESPONSE')
xlabel('Frequency (yr^{-1})')
ylabel('Amplitude')
text(0.25,0.70,['90% Response: ',num2str(p90),' yr']);
text(0.25,0.65,['50% Response: ',num2str(p50),' yr']);
text(0.25,0.60,['10% Response: ',num2str(p10),' yr']);


