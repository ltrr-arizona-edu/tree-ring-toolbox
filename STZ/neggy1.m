% neggy1.m  -- try negative exponential

dummy=3;
notan=NaN;
bb = [.50 .50];

figure(1)
t=(1:500)';
k=1.0
a= 0.8;
b= -0.04;
y=k + a*exp(t*b)  + 0.1*randn(500,1);

W=[t y];
xlab=' ';
ylab=' ';
nm = 'this';
x0=0.1;


%[x] = leastsq('negxpk1',x0,[],[],W);
%[f,g,c]=negxpk1(x,W);

yrv =t+1200;
yrvn=t+1200;
xvn = y;
xvn(120:140)=dummy(ones(21,1),:);

%xvn(1:21)=notan(ones(21,1),:);
%yrvn(1:21)=notan(ones(21,1),:);

xvn(120:140)=notan(ones(21,1),:);
yrvn(120:140)=notan(ones(21,1),:);

[cvx,k,a,b] = cfnegx(yrv,yrvn,xvn)
