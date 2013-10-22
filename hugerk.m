function f=hugerk(x,w)

n=length(w);
for t=1:n,
	g(t)= (x(1)*t .^x(2)) * exp(-t*x(3))  + x(4);
	f(t)=w(t)-g(t);
end
	