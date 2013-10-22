function y=fun1(x,p)
f=x;
A=6*(cos(2*pi*f)-1).^2;
B=(cos(2*pi*f) + 2);

y = (A/B - p).^2;
