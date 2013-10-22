% endyrs.m

% forms two vectors, (go and stop), containing the first and
% last valid years of each column beyond column 1 of a loaded
% input array A.

% saves these vectors in the file limyrs.mat

% assumes that invalid data is coded as zeros before and after any
% valid data for series.  This can be changed easily.

% used on MKH's sequoia array to prepare for pdgm1.m
% assumes the time-series array A has been pre-loaded


ago=input('First row of A corresponds to what year?  ');

for n=2:13
	I=find(A(:,n)>0);
	go(n-1)=min(I)+ ago -1;
	stop(n-1)=max(I) + ago -1;
	clear I
end

save limyrs go stop
