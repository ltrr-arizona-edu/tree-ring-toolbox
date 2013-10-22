function tsmcorei(IX,EX,yrs,nms)
% Make time series matrix of core indices


missv = 9990;
s1 = input('Standard (S) or Residual (R) Indices? ','s')

if s=='S',
	X = IX;
	T = yrs;
elsif s =='R',
	X = EX;
	T=[yrs(:,1)+arord  yrs(:,2)  yrs(3)+arord];
end




