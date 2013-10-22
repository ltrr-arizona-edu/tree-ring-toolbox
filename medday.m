function y2=medday(y1,x1,x2)
% median ratio estimation of missing ppt data; called by dayest.m
%
% D. Meko 11-8-93
%
%
% y1 --- cv of ppt values for key (predictand) station
% x1 --- cv of contemporaneous values for predictor station
% x2 --- cv of values at predictor station for same times that data
%		missing at key station
% 
% y2 --- estimated cv for key station corresponding to values in x2
%
%
%**************


% Assumes no zero or negative values in y1, x1 or x2 


mx2=length(x2);
y2=zeros(mx2,1);

% check for missing or zero values
L1=[any(y1==0)  any(x1==0)  any(x2==0)];
if any(L1);
	error('zero ppt passed to medday.m in x1,y1 or x2');
end
L2=[any(y1<0)  any(x1<0)  any(x2<0)];  
if any(L2);
	error('negative ppt passed to medday in x1, y1 or x2');
end

% Compute median ratio
z= x1 ./ y1;  % ratio series
md = median(z);  % median ratio

% Multiply median ratio times predictor series to get estimates
y2 = md * x2;
end
