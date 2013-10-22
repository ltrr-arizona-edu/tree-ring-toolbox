M=15;  % lag

for i=1:3
	x=B(:,i+1);
	x=x-mean(x);
	cv=covf(x,M);  % autocovariance function
	av=cv ./ cv(1);

	y(:,i)=av';
end
