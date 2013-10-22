function b=lowpass1
%
% FIR filter design using matlab function fir1
%
% Meko 7-24-96

clc
disp('Number of weights must be odd');
nw = input('Number of weights in filter');

fhalf=input('Frequency to try for amplitude of 0.5:  ');


% Convert frequency into desired units
whalf=fhalf*2;

% Get filter 
b = fir1(nw-1,whalf);


% Plot filter weights in figure 1
figure(1)
plot(b);
title('Filter Weights')




% Compute frequency response
% h is the complex response
% w is the frequency in radians (pi radians is 0.5 in f units)
[h,w]=freqz(b,1,100);


% Convert frequency vector w to  [0 .5] units
f = w /(2*pi);


% Get amplitude of frequency response
a = abs(h);

% Getphase
ph=angle(h);


% plot amplitude and phase
figure(2)
%subplot(2,1,1);
plot(f,a);
title('Amplitude of Frequency Response');
%subplot(2,1,2);
%plot(f,ph);
%title('Phase')

