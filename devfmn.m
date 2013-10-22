function [P,m] = devfmn(X)

%  Transform a data array X into deviations from long-term mean
%  Return the col means in m.

% Used on temperature data to convert to departures

[m1,n1]=size(X);
m=mean(X);
M=m(ones(m1,1),:);
P=X - M;
