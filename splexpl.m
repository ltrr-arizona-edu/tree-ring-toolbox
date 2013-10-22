% splexpl.m

% Build time series, y, length 530 yr, composed of a 53-yr sin
% Compute spline-smoothed version, s,  using p=1e-4 
% Compute "index" series z by ratio y/s

% Compare variances:


t=1:530;
w=2*pi*t/53;
y=sin(w);
y=y+1;    % to make y non-negative

s = csaps(t,y,1e-4,t);

z=y ./ s;


plot(t,y,t,s,t,z)
pause


ss=[std(y) std(s)  std(z)];  % standard deviations
vv=ss .^2;
disp(vv)
