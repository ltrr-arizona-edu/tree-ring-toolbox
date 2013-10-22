function [NSTNS,NS,yr3]= regsub1(ns1,nreg,F,J2,IR1,yr2)
%
%
% Sample size (number of stations) in regional climate calculations
% Dedicated function regcli2.m
%
% Meko 2-8-97
%
% Input args defined in regcli2.m
% Output args
%  NSTNS number of station matrix for overall study (all regions)
%  NS  likewise, but for each region. Cols 1-12 for reg 1, etc
%  yr3 cv of years for NSTNS and NS

nyrs2=length(yr2);

NSTNS=repmat(NaN,nyrs2,12);
if ~isempty(nreg);
   NS=repmat(NaN,nyrs2,(12*nreg));
else;
   NS=[];
end;

% Logical matrix, 1 if data this station/month/year, 0 otherwise
L1=~isnan(F);

%**************** FIRST TALLY IS IRRESPECTIVE OF REGION *************


for n=1:12; % loop over months
	jcols=J2(:,n); % index to cols of F (or L1) for this month of yr
	LL1=L1(:,jcols); % mtx of elements of L1 for this month/yr
	n1=(nansum(LL1'))';
	NSTNS(:,n)=n1;
end


%**** NEXT TALLY GIVE CHANGE IN SAMPLE SIZE BY REGION ***********

if ~isempty(nreg);
   for n=1:12; % loop over months  of year
      % pull this month's data for all ns1 stns
      jcols=J2(:,n);
      LL1=L1(:,jcols);
      % loop over regions
      for k=1:nreg
         % calc pointer for result cv to go into NS
         j1=n+(k-1)* 12;
         % pull elements of jcols for this region's stations
         ir1=IR1(:,k); % station number in this region -- zero padded cv
         ir1(ir1==0)=[]; % get rid of zero elements
         LL2=LL1(:,ir1); % logical matrix, this region's stations
         n2=(sumnan(LL2'))'; % count for this month/region
         NS(:,j1)=n2; % store result
      end
   end
end;

%********** TRIM AL-NAN ROWS FROM NS, NSTNS

yr3=yr2; % initialize year vector for trimmed mtx
L1=NSTNS==0;
L2=(all(L1'))';
if sum(L2)>0;
   NSTNS(L2,:)=[];
   if ~isempty(nreg);
      NS(L2,:)=[];
   end;
   yr3(L2)=[];
end

% Just in case, check that yr3 continuous
d1=diff(yr3);
if  ~all(d1==1);
   yr3
   error('yr3 not continous');
end;

