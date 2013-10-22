function  hist2(mv,oc,ov)
subplot(211)
L=min([min(oc) min(ov)]);
M=max([max(oc) max(ov)]);
d=(M-L)/10;
x=linspace((L-d):(M+d),10);

hist(oc,x)
hist(ov,x)
pause
subplot(111)
