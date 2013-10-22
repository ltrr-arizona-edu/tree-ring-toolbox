function f = fitfun3(lam,Data)
%FITFUN3  Used  to return errors in fitting data to
%	  Hugershoff function with nonzero asymptote. Errors are 
%	between the data and the values computed by the current
%	function of lam.  FITFUN3 assumes a function of the form
%
%	  y = c(1)+ c(2)* (t .^lam(1)) .* exp(-lam(2)*t) 
%
%	with 2 linear parameters and 2 nonlinear parameters.

n=length(Data);
t = (1:n)';
y=Data;
A = zeros(n,2);
A(:,1) = ones(n,1);
A(:,2) = (t .^lam(1)) .* exp(-lam(2)*t);

c = A\y;
z = A*c;
f = z-y;

% Statements to plot progress of fitting:
plot(t,z,t,y,'o')
xt = max(t)/3;
yt = max(y)/2;
text(xt,1.4*yt,['lambda = ' num2str(lam(1)) '  ' num2str(lam(2))])
text(xt,1.2*yt,['(c = ' num2str(c(1)) '  ' num2str(c(2)),')'])
text(xt,1.0*yt,['err norm = ' sprintf('%g',norm(f))])
