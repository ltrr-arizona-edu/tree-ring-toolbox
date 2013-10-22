
% fftcheck.dat

% load pdf11
xx=dtrend(xa,0);
x=zeros(256,1);
x(1:length(xa))=xx;

y=fft(x,256)

