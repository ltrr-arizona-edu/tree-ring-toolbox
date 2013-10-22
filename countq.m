function [q,k]= countq(x,p)

% count number of years outside quantiles in subperiod

%***************   USER-WRITTEN FUNCTIONS NEEDED
%
% quantile.m  --- to compute quantiles
% subyears.m  --- to compute beg, end years of sub periods
% 


y1on = input('Start year for entire time series:  ');
y1off = input('End year for entire time series:  ');

y2on = input('Start year for modern period:  ');
y2off = input('End year for modern period:  ');

yr=(y1on:y1off)';

k= subyears([y1on y1off], [y2on  y2off]);  % compute subperiods
q=quantile(x,p);

[mk,nk]=size(k);

S=zeros(mk,1);

for i=1:mk;
	L(:) = yr>= k(i,1)  & yr <= k(i,2);
	if (p<=0.5);  %  dry side
		Lq= x<=q;
	else
		Lq= x>=q;
	end
	L2=L & Lq;
	S(i) = sum(L2);
end

k=[k S];


plot(yr,x,[yr(1)  yr(length(yr))], [q q]);
pause
