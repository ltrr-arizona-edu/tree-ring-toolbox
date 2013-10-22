% CODE FOR GENERATING TS PLOTS OF 3 TREE-RING SERIES

V=[1700 2000 0 2];
axis(V)
subplot(221)
yr=(1700:1979)';
plot(yr,ca,'-');
plot(yr,ke,'-');
plot(yr,ga,'-');


% CODE FOR HISTOGRAMS

clg
subplot(111)

subplot(221)

hist(ca);
text(.6,50,'Cass Lake, PIRE')

hist(ke)
text(.6,50,'Ketchum, PIFL')

hist(ga)
text(.6,50,'Gallinas, PIED')



%  CODE FOR GETTING PLOTS OF ACFS

M=15;  % lag

for i=1:3
	x=B(:,i+1);
	x=x-mean(x);
	cv=covf(x,M);  % autocovariance function
	av=cv ./ cv(1);

	y(:,i)=av';
end
