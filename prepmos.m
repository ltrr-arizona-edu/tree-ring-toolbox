function [A,mlabs]=prepmos(P,T,yrs,endmo,nmos)
% prepmos: prepare lagged monthly matrix of precipitation and temperature
% [A,mlabs]=prepmos(P,T,yrs,endmo,nmos);
% Last revised 7-13-99
%
% Form up-to-18-months lagged matrix of monthly ppt and temperature.  Called by 
% other functions (e.g., respfun1.m) that deal with the response of tree rings 
% to monthly and seasonal climate.
%
%*** INPUT ARGS
%
% P (? x 13)r  monthly precip series, year as col 1.  No missing data.
% T (? x 13)r  monthly temp   series, year as col 1.  No missing data.
% yrs (1 x 2)i  beginning and ending year of desired output lagged array.
% endmo (1 x 1)  ending month for climate year:  1=jan  12=dec
% nmos (1 x 1)  number of months (<=18) in climate year
%
%
%*** OUTPUT ARGS
%
% A  the lagged ppt and temperature array
%   number of rows = yrs(2) - yrs(1) + 1
%   number of cols =  nmos*2 + 1
%   First cols are ppt, then temp (see mlabs)
% mlabs (nmos*2 x 1)  labels for columns of A (not incl year)
%   Convention: P00J -- precip jan of current year
%               P00F -- precip feb of current year
%               P-1J -- precip jan of year t-1
%               T-1J -- temp   jan of year t-1
%                  etc
%
%*** REFERENCES -- none
% 
%*** UW FUNCTIONS CALLED -- none
%*** TOOLBOXES NEEDED -- none
%
%*** NOTES
%
% The years in yrs are defined relative to the ending month of the period as specified 
% by endmo.  For example, if yrs is [1900 1990], endmo is 8, and nmos is 18, the 
% period begins with the 18-months April 1899-September 1900, and ends with the 18
% months April 1989-September 1990.

[mp,np]=size(P);    % number row, cols in input precip array
[mt,nt]=size(T);    % same for monthly mean temp.

%*** Calc number of side-by-side blocks of data needed for lagged array
%  If nblks = 0 monthly grouping can be satisfied with data from current
%  yr.  If nblks=1, need one previous year.  If nblks=2, need 2 previous

if nmos-endmo >= 13
	nblks=2;
elseif nmos-endmo < 1
	nblks=0;
else
	nblks=1;
end


%**************  CHECK THAT P, T HAVE REQUIRED YEAR COVERAGE TO COVER
%     TARGET OUTPUT PERIOD -- INCLUDING LAGS

L3 = [P(1,1)>yrs(1)-nblks   P(mp,1)<yrs(2)  T(1,1)>yrs(1)-nblks ...
      T(mt,1)<yrs(2)];
if any (L3)
	error('INSUFFICIENT DATA COVERAGE IN P OR T FOR SPECIFIED YRS');
end



%************  CHECK THAT NO YEARS ARE SKIPPED IN THE P, T ARRAYS

nyrs = yrs(2)-yrs(1) +1;  % Desired output number of years
L1= (P(:,1) >= yrs(1)  & P(:,1) <= yrs(2)); % 0/1 pointer to years in P
L2= (T(:,1) >= yrs(1)  & T(:,1) <= yrs(2)); %... and in T
if any([sum(L1)-nyrs   sum(L2)-nyrs]);
	error('DATA INPUT IN EITHER P OR T SKIPS A YEAR ');
end

%*********   SET POSSIBLE LABELS FOR MONTHS

namp1=['P-2J';'P-2F';'P-2M';'P-2A';'P-2M';'P-2J'];
namp2=['P-2J';'P-2A';'P-2S';'P-2O';'P-2N';'P-2D';'P-1J'];
namp3=['P-1F';'P-1M';'P-1A';'P-1M';'P-1J';'P-1J';'P-1A';'P-1S'];
namp4=['P-1O';'P-1N';'P-1D';'P00J';'P00F';'P00M';'P00A';'P00M';'P00J'];
namp5=['P00J';'P00A';'P00S';'P00O';'P00N';'P00D'];
namp=[namp1;namp2;namp3;namp4;namp5];


namt1=['T-2J';'T-2F';'T-2M';'T-2A';'T-2M';'T-2J';'T-2J';'T-2A';'T-2S'];
namt2=['T-2O';'T-2N';'T-2D';'T-1J';'T-1F';'T-1M';'T-1A';'T-1M';'T-1J'];
namt3=['T-1J';'T-1A';'T-1S';'T-1O';'T-1N';'T-1D';'T00J';'T00F';'T00M'];
namt4=['T00A';'T00M';'T00J';'T00J';'T00A';'T00S';'T00O';'T00N';'T00D'];
namt=[namt1;namt2;namt3;namt4];

ne=endmo+24;
nb=ne-nmos+1;
mlabs=[namp(nb:ne,:) ; namt(nb:ne,:)];

%*******  

P2=sideby(P,L1,endmo,nblks,nmos);
T2=sideby(T,L2,endmo,nblks,nmos);

A =[P(L1,1)   P2   T2];


%*** SUBFUNCTIONS

function P2=sideby(P,L1,endmo,nblks,nmos)
%
% sub-function of prepmos.m
%
%****   INPUT ARGS
%
%  P (? X 13)   array of montly ppt or temp
%  L1 (? x 1)  0-1 pointer to specific years of P
%  endmo (1 x 1)  ending month of year for climate analysis
%  nblks (1 x 1)  number of blocks to stack side by side
%     nblks is computed within calling pgm  prepmos.m
%
%*** OUTPUT ARGS
%
%  P2 (? x ?)  output lagged array

L1L = logical([L1(2:length(L1)); 0]);
L1LL= logical([L1(3:length(L1)); 0; 0]);

P1=P(L1,2:endmo+1);
if nblks==1
	P2 = [P(L1L,2:13)  P1];
	[mp2,np2]=size(P2);
	ncout = np2-nmos;
	P2(:,1:ncout)=[];
elseif nblks==2
	P2=[P(L1LL,2:13)  P(L1L,2:13)  P1];
	[mp2,np2]=size(P2);
	ncout=np2-nmos;
	P2(:,1:ncout)=[];
else
	P2=P1;
	[mp2,np2]=size(P2);
	ncout=np2-nmos;
	P2(:,1:ncout)=[];
end

