% neggy.m  -- try negative exponential


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


[x] = leastsq('negexpk',x0,[],[],W,nm,xlab,ylab);
pause

%figure(2)
%z=filtfilt(bb,1,y);
%W=[t,z];
%nm='that';

%[xx] = leastsq('negexpk',x0,[],[],W,nm,xlab,ylab);


[f,g,c]=negexpk(x,W,nm,xlab,ylab)
