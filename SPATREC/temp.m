function e=pltrec3(yr4,modnum,EV,RE,NV)
% e=pltrec3(yr4,modnum,EV,RE,NV)

a=NaN;

% Compute time series of EV,RE,n1,n2,n3
for n = 1:length(modnum);
	zev(n) = EV(modnum(n));
	zre(n) = RE(modnum(n));
	zn1(n) = NV(modnum(n),3); % available trees
	zn2(n) = NV(modnum(n),4); % potential pc-amp predictors
	zn3(n) = NV(modnum(n),5); % final pc-amp predictors
end

subplot(2,1,1)
h1=plot(yr4,EV,yr4,RE)

v1=axis;
ylabel(' ')
xlabel('Year')
title('EXPLAINED VARIANCE AND REDUCTION OF ERROR')


subplot(2,1,2)
plot(yr4,zn1,yr4,zn2,yr4,zn3);
axis(v1)
ylabel('Number of trees')
legend('No. of Trees','No. Potential Predictors','Number Predictors')


