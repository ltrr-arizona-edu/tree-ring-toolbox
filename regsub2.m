function [NS,yr3]= regsub1(ns1,F,J2,yr2)
%
%
% Sample size (number of stations) in target point estimation each year
% Dedicated function regcli3.m
%
% Meko 11-20-97
%
% Input args defined in regcli3.m
% Output args
%  NS  valid sample size for each year/month
%  yr3 cv of years for NS

nyrs2=length(yr2);

a=NaN;

NS=a(ones(nyrs2,1),ones(12,1));

% Logical matrix, 1 if data this station/month/year, 0 otherwise
L1=~isnan(F);

%**************** FIRST TALLY IS IRRESPECTIVE OF REGION *************


for n=1:12; % loop over months
	jcols=J2(:,n); % index to cols of F (or L1) for this month of yr
	LL1=L1(:,jcols); % mtx of elements of L1 for this month/yr
	n1=(sumnan(LL1'))';
	NS(:,n)=n1;
end


%********** TRIM ALL-NAN ROWS FROM NS, NSTNS

yr3=yr2; % initialize year vector for trimmed mtx
L1=NS==0;
L2=(all(L1'))';
if sum(L2)>0;
	NS(L2,:)=[];
	yr3(L2)=[];
end


% Just in case, check that yr3 continuous
d1=diff(yr3);
if  ~all(d1==1);
	yr3
	error('yr3 not continous')
end

