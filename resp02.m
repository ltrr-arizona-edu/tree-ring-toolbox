function [EV,SS]=resp02(miscfl)
% resp02:  Output-error (OE) modeling of tree-ring response to monthly climate
% CALL: [EV,SS]=resp02(miscfl);
%
% Meko  10-28-97
%
% Application:  Screening of a tree-ring index against monthly
% precipitation and temperature data to find out the best seasonal 
% grouping of months for a climate signal.  Produces colormaps summarizing
% strength of signal and direction (+ or -) of relationship. Handles
% distributed lag response by fitting OE models of different orders: (1,0),
% (2,0) and (1,1).  Signal strength measured by 1 - [var(errors)/var(index)].
%
%
%*******************  INPUT ********************************
%
% miscfl (1 x ?)s name of .mat file with the following variables
%
%		yrgo (1 x 1)i -- start year of analysis period
%		yrsp(1 x 1)i -- end year of analysis period
%		endmo(1 x 1)i -- ending month of growth year
%		xmiss (1 x 1)r-- missing value numeric code
%		kopt(1x1)i -- options: reserved; for now, just set as 1
%    Afile(1 x ?)s  name of .mat file with A data (usually pcp) -- assumed in a 
%       matrix X,Y,or Z
%    Bfile(1 x ?)s  name of .mat file with B data (usually tmp) -- assumed in a
%       matrix X,Y, or Z
%    Cfile (1 x ?)s name of .mat file with tree-ring index -- assumed in a 
%       matrix X, Y or Z
%		varnma (1 x ?)s code for A variable (e.g., 'PCP') -- max of 6 chars
%		varnmb (1 x ?)s code for B variable (e.g., 'TMP') -- max of 6 chars
%		sitecode (1 x ?)s site code (e.g., 'KCSSTD') -- max of 6 chars
%
%***************   OUTPUT  ARGS  *******************************************
%
%
% EV{} cell variable with explained variance proportion for models
%
%   EV{1} (12 x 12 x 3)r results for modeling of tree vs the 'A' variable,
%		typically pcp. A given row represents a specific ending month of a monthly
%		grouping. A given column represents a specific number of months in the grouping. 
%		Dimension 3 represents the three different orders of OE models:
%		1 is OE(1,0);  2 is OE(2,0);  3 is OE(1,1)
%
%	  EV{2} same as EV{1}, except for tree vs the 'B' variable, typically
%			temperature
%
% SS{} coded sign and importance of response of growth to climate.
%   SS{1} (12 x 12 x3)i  ppt, 
%			15==sign of relationship (sum of first few irf weights) negative, 
%           all B and F model parameters at least 3 std devs from zero, and
%				OE model explains at least 5% of tree-ring variance
%			45==Like '15' case above, except sign of relationship positive
%        30==Relationship weak:  either explained variance less than 5%, or
%           all model parameters not at least 3 std devs from zero
%		SS{2} ... like SS{1}, but for tmp
%
%
%***************** OUTPUT, GRAPHICAL
%
% Pairs of figure windows with color maps are produced.  First pair (windows 1 and 2)
% contain results for OE(1,0) model on PCP. Window 1 color maps the signal strength, 
% as measured by 1 - [var(errors)/var(x)], where x is the output tree-ring index
% Window 2 maps the sign of the relationship between PCP and the
% tree-ring index for those monthly groupings with a 'strong' signal
% Figs 3,4: like above, except OE(2,0) model
% Figs 5,6: like above, except OE(1,1) model
% Figure windows 7-12:  like windows 1-6, but for TMP instead of PCP
%
%
%*************   PRELIMINARY STEPS  ****************************************
%
% Get the pcp and tmp data into .mat files with year in col 1 and Jan-Dec
% data in cols 2-13.  Store the data in matrices named X, Y or Z.
% Pcp and tmp matrices need not have identical year coverage, but both
% must have non-missing data for the specified analysis period yrgo:yrsp, and
% for the 2 years preceding yrgo
% Missing values must be NaN or a missing value code specified by 
% xmiss
%
% Get the tree-ring standard chronology into a matrix, with year in
% col 1 and the index in col 2.  Store in matrix named X, Y and Z.  All
% data for analysis period yrgo:yrsp must be valid
%
% Prepare the miscellaneous input .mat file that has the file names of
% the .mat files with pcp, tmp and tree-ring data, as well as other
% model specifications
%
%**************** NOTES **********************************
%
% In a layer of EV{1},EV{2},SS{1} or SS{2}, element (i,j) holds statistic for
% the j-month period ending in 'relative' month i -- that is
% relative to the tree-year end month specified by endmo.
% For example, if endmo=8 (August)
% row 1 holds results for all monthly groupings ending in August of year t;
% row 2 is for groupings ending in July of year t ..., 
% and row 12 is for groupings ending in September of  year t-1.
%
% resp02.m written in terms of precipitation (PCP) and temperature (TMP) as
% the primary and secondary climatic variables.  But these could optionally
% be other variables (e.g., soil moisture and wind speed)
%
%******************** END OF OPENING COMMENTS ******************************
clc



%***********************  LOAD INPUT FILES



%----------  MISCELLANEOUS INPUT


eval(['load ' miscfl ';']);
Ltemp = [exist('yrgo')==1  exist('yrsp')==1  exist('endmo')==1 ...
      exist('kopt')==1 exist('xmiss')==1  exist('Afile')==1 ...
      exist('Bfile')==1  exist('Cfile')==1 exist('varnma')==1 ...
		exist('varnmb')==1  exist('sitecode')==1];
if ~all(Ltemp);
	error('miscfl does not have all required input variables');
end


%--------- A FILE (USUALLY PCP)

eval(['load ' Afile ';']);
if exist('X')==1;
	A=X;
	clear Y Z;
elseif exist('Y')==1;
	A=Y;
	clear X Z;
elseif exist('Z')==1;
	A=Z;
	clear X Y;
else
	error(' Afile does not contain X, Y, or Z');
end
clear X Y Z


%--------- B FILE (USUALLY TMP)

eval(['load ' Bfile ';']);
if exist('X')==1;
	B=X;
	clear Y Z;
elseif exist('Y')==1;
	B=Y;
	clear X Z;
elseif exist('Z')==1;
	B=Z;
	clear X Y;
else
	error(' Bfile does not contain X, Y, or Z');
end
clear X Y Z


%--------- C FILE (USUALLY tree-ring index)

eval(['load ' Cfile ';']);
if exist('X')==1;
	C=X;
	clear Y Z;
elseif exist('Y')==1;
	C=Y;
	clear X Z;
elseif exist('Z')==1;
	C=Z;
	clear X Y;
else
	error(' Cfile does not contain X, Y, or Z');
end
clear X Y Z


%********************  CHECK YEARS AND BUILD YEAR POINTERS

yrC = C(:,1);   % col vector of years
yrsC=[min(yrC)  max(yrC)]; % start and end year
yrA = A(:,1);
yrsA=[min(yrA)  max(yrA)];
yrB = B(:,1);
yrsB=[min(yrB)  max(yrB)];
yr1 = (yrgo:yrsp)';  % ........ for desired analysis period
yrs1 = [min(yr1) max(yr1)];


% Check starting years, keeping in mind that will need 2 leading years of
% data on the A and B matrices
if yrsA(1)>yrgo-2 | yrsB(1)>yrgo-2 | yrsC(1)>yrgo;
	error('A, B, or C has too late a starting year for yrgo');
end
% Check ending years
if yrsA(2)<yrsp | yrsB(2)<yrsp | yrsC(2)<yrsp;
	error('A, B, or C has too early an ending year for yrsp');
end


%-------------- MENU FOR WHETHER TO DO BOTH A AND B VARIABLE, OR JUST A
kmods=menu('MODEL JUST ONE VARIABLE OR BOTH?',...
   ['Both ' varnma ' and ' varnmb],...
   [varnma ' only'],...
   [varnmb ' only']);
switch kmods;
case 1
   kvar = 1; % model both
case 2
   kvar=2;  % Model just A 
case 3
   kvar=3; % Model just B
otherwise
   error('Illegal entry for kmods');
end



%-------------- DIALOG FOR AUTOMATIC PRODUCTION OF .PS FILES
kpost = menu('Automatically build a .PS file for each figure window?',...
   'Yes',...
   'No');
if kpost==1;
   kpost='Yes';
else;
   kpost='No';
end

%kpost=questdlg('Automatic Building of .PS file for each figure window?');
% kpost is now 'Yes', 'No', or 'Cancel'



%------------------ YEAR POINTERS

LAk = yrA>=yrgo-2 & yrA<=yrsp;  % keeper data for A -- should have no missing
LBk = yrB>=yrgo-2 & yrB<=yrsp;  % keeper data for B
LC = yrC>=yrgo & yrC<=yrsp; % to analysis period in C
LA = yrA>=yrgo & yrA<=yrsp; % to analysis period in A
LA1 = yrA>=yrgo-1 & yrA<=yrsp-1 ; % to lag minus 1 in A
LA2 = yrA>=yrgo-2 & yrA<=yrsp-2;  % to lag minus 2 in A
LB = yrB>=yrgo & yrB<=yrsp; % to analysis period in B
LB1 = yrB>=yrgo-1 & yrB<=yrsp-1 ; % to lag minus 1 in B
LB2 = yrB>=yrgo-2 & yrB<=yrsp-2;  % to lag minus 2 in B

% Check for missing data 
if any(any(A(LAk,:)==xmiss));
	error('A has missing data in analysis period or 2 leading years');
end
if any(any(B(LBk,:)==xmiss));
	error('B has missing data in analysis period or 2 leading years');
end
if any(any(C(LC,:)==xmiss));
	error ('C has missing data in analysis period');
end
if any(any(isnan(A(LAk,:))));
	error('A has missing data in analysis period or 2 leading years');
end
if any(any(isnan(B(LBk,:))));
	error('B has missing data in analysis period or 2 leading years');
end
if any(any(isnan(C(LC,:))));
	error ('C has missing data in analysis period');
end


%*******************  BUILD LAGGED MATRICES C AND D BASED ON A AND B

% Year is in col 1, lag minus 2 data in next  12 cols, minus-1 in next,
% current years data in last 12 cols
D = [yr1 A(LA2,2:13)  A(LA1,2:13)  A(LA,2:13)]; % PCP matrix
E = [yr1 B(LB2,2:13)  B(LB1,2:13)  B(LB,2:13)]; % TMP matrix



%****************  Pull off analysis period years of tree-ring series
C2=C(LC,2);




%*************** OUTPUT ERROR MODELING


t1= 'The modeling has begun.  Typically takes 2-3 minutes.  Resp02.m ';
t1=str2mat(t1,'first fit all the models to the PCP variable, then to the TMP '); 
t1=str2mat(t1,'variable.  Twelve figure windows will then be created.  Windows ');
t1=str2mat(t1,'1-6 map results for PCP and windows 7-12 for TMP.  For PCP, Windows');
t1=str2mat(t1,'1,3,5 map the strength of signal in terms of tree-ring variance ');
t1=str2mat(t1,'accounted for by the OE(1,0), OE(2,0) and OE(1,1) models, and windows');
t1=str2mat(t1,'2,4,6 give corresponding maps of sign of the relationship, colored for');
t1=str2mat(t1,'those variables whose OE models have at least 1 parameter at least');
t1=str2mat(t1,'3 standard deviations from zero.  Windows 7,9,11 and 8,10,12 similarly');
t1=str2mat(t1,'map results for the TMP variable.');
disp(t1);
disp('  ');
disp('Press any key to continue with the modeling')
pause


%-----------------  tree vs pcp
if kvar==1 | kvar==2;
   disp('Modeling tree rings against pcp ...');
   [F,G]=respoe(C2,D,endmo,1,1); %  subfunction that does the modeling
   EV{1}=F;
   SS{1}=G;
end


%-----------------  tree vs tmp
if kvar==1 | kvar==3;
   disp('Modeling tree-rings against tmp...');
   [F,G]=respoe(C2,E,endmo,2,1);
   EV{2}=F;
   SS{2}=G;
end



%******************* FIRST-ORDER AUTOCORRELATION COMPUTATION (NOT YET USED)

% First-order autoc of pcp
if kvar==1 | kvar==2;
   %disp('Modeling autocorrelation of pcp');
   F=respoe(C2,D,endmo,1,2);
   EV{3}=F;
end


% First-order autoc of tmp
if kvar==1 | kvar==3;
 %disp('Modeling autocorrelation of tmp');
 F=respoe(C2,E,endmo,2,2);
 EV{4}=F;
end



%****************  MAP RESULTS

disp('BUILDING FIGURE WINDOWS WITH GRAPHS OF RESULTS');
    
% Build a yellow, white green color map for the 'sign of signal' windows.
% Red, white or green will be the colors, depending on sign of relationship and
% strength.
wgy = [repmat([1 0 0],20,1); repmat([1 1 1],20,1) ; repmat([0 1 0],24,1)];
    
    
txtin{2}='Number of Months in Window';
txtin{3}='Ending Month of Window';
txtin{4}='jet';
txtyr = sprintf('%4.0f-%4.0f',yrgo,yrsp);
datin{1}=endmo;
datin{3}=1;

% Set color limits for signal strength to max and min of any variance explained in any of the 
% 144 x 3 possible models (12 ending months x 12 numbers of months by 3 
% model orders)
cmaxA = max(max(max(EV{1})));
cmaxB = max(max(max(EV{2})));
cminA = min(min(min(EV{1})));
cminB = min(min(min(EV{2})));

if kvar==1; % set color limits from results on both types of variables
   cmin = min(cminA,cminB);
   cmax = max(cmaxA,cmaxB);
elseif kvar==2; % pcp only
   cmin = cminA;
   cmax = cmaxA;
elseif kvar==3; % tmp only
   cmin = cminB;
   cmax = cmaxB;
end
clim = [cmin cmax];
datin{4}=clim;


%------- PCP MAPPING
if kvar==1 | kvar==2;
   G=SS{1};
	Y = EV{1};
   txt1 = varnma;
   
   %--------- OE(1,0)
   Y1 = Y(:,:,1);
   G1  = G(:,:,1);
	txt2='OE(1,0)';
   nfig1=1; nfig2=2; % figure windows
   
   txtin{1}=['Signal Strength -- ' txt1 '--' txt2 ' Model -- ' sitecode];
	txtin{1}=[txtin{1} ' (' txtyr ')'];
   datin{2}=flipud(Y1);
   figure(nfig1);
   colmap01(txtin,datin);
   if strcmp(kpost,'Yes'); % postscript file
      eval(['print -dpsc ' 'wind' int2str(nfig1) ';']);
   end

   datina = datin;
   txtina=txtin;
   txtina{4}=wgy;
   txtina{1}=['Sign of Response -- ' txt1 '--' txt2 ' Model -- ' sitecode];
	txtina{1}=[txtina{1} ' (' txtyr ')'];
   datina{2}=flipud(G1);
   datina{4}=[1 64];
   figure(nfig2);
   colmap02(txtina,datina);
   if strcmp(kpost,'Yes'); % postscript file
      eval(['print -dpsc ' 'wind' int2str(nfig2) ';']);
   end
   
   %------------ OE(2,0)
   
   Y1= Y(:,:,2);
   G1= G(:,:,2);
	txt2='OE(2,0)';
   nfig1=3;  nfig2=4;
   
   txtin{1}=['Signal Strength -- ' txt1 '--' txt2 ' Model -- ' sitecode];
	txtin{1}=[txtin{1} ' (' txtyr ')'];
   datin{2}=flipud(Y1);
   figure(nfig1);
   colmap01(txtin,datin);
   if strcmp(kpost,'Yes'); % postscript file
      eval(['print -dpsc ' 'wind' int2str(nfig1) ';']);
   end

   datina = datin;
   txtina=txtin;
   txtina{4}=wgy;
   txtina{1}=['Sign of Response -- ' txt1 '--' txt2 ' Model -- ' sitecode];
	txtina{1}=[txtina{1} ' (' txtyr ')'];
   datina{2}=flipud(G1);
   datina{4}=[1 64];
   figure(nfig2);
   colmap02(txtina,datina);
   if strcmp(kpost,'Yes'); % postscript file
      eval(['print -dpsc ' 'wind' int2str(nfig2) ';']);
   end
   
   
   %-----------  OE(1,1) model
   Y1=Y(:,:,3);
   G1=G(:,:,3);
	txt2='OE(1,1)';
	nfig1=5; nfig2=6;
   
   txtin{1}=['Signal Strength -- ' txt1 '--' txt2 ' Model -- ' sitecode];
	txtin{1}=[txtin{1} ' (' txtyr ')'];
   datin{2}=flipud(Y1);
   figure(nfig1);
   colmap01(txtin,datin);
   if strcmp(kpost,'Yes'); % postscript file
      eval(['print -dpsc ' 'wind' int2str(nfig1) ';']);
   end

   datina = datin;
   txtina=txtin;
   txtina{4}=wgy;
   txtina{1}=['Sign of Response -- ' txt1 '--' txt2 ' Model -- ' sitecode];
	txtina{1}=[txtina{1} ' (' txtyr ')'];
   datina{2}=flipud(G1);
   datina{4}=[1 64];
   figure(nfig2);
   colmap02(txtina,datina);
   if strcmp(kpost,'Yes'); % postscript file
      eval(['print -dpsc ' 'wind' int2str(nfig2) ';']);
   end
end

if kvar==1 | kvar==3;  % variable B (e.g., temperature)
   
   figure(2)
   Y = EV{2};
   G = SS{2};
   txt1 = varnmb;
      
   %--------  OE(1,0) 
   Y1 = Y(:,:,1);
   G1=G(:,:,1);
	txt2='OE(1,0)';
   nfig1=7; nfig2=8;
   
   txtin{1}=['Signal Strength -- ' txt1 '--' txt2 ' Model -- ' sitecode];
	txtin{1}=[txtin{1} ' (' txtyr ')'];
   datin{2}=flipud(Y1);
   figure(nfig1);
   colmap01(txtin,datin);
   if strcmp(kpost,'Yes'); % postscript file
      eval(['print -dpsc ' 'wind' int2str(nfig1) ';']);
   end

   datina = datin;
   txtina=txtin;
   txtina{4}=wgy;
   txtina{1}=['Sign of Response -- ' txt1 '--' txt2 ' Model -- ' sitecode];
	txtina{1}=[txtina{1} ' (' txtyr ')'];
   datina{2}=flipud(G1);
   datina{4}=[1 64];
   figure(nfig2);
   colmap02(txtina,datina);
   if strcmp(kpost,'Yes'); % postscript file
      eval(['print -dpsc ' 'wind' int2str(nfig2) ';']);
   end

   %---------- OE (2,0)  model
   
   Y1= Y(:,:,2);
   G1=G(:,:,2);
	txt2='OE(2,0)';
   nfig1=9; nfig2=10;
   
   txtin{1}=['Signal Strength -- ' txt1 '--' txt2 ' Model -- ' sitecode];
	txtin{1}=[txtin{1} ' (' txtyr ')'];
   datin{2}=flipud(Y1);
   figure(nfig1);
   colmap01(txtin,datin);
   if strcmp(kpost,'Yes'); % postscript file
      eval(['print -dpsc ' 'wind' int2str(nfig1) ';']);
   end

   datina = datin;
   txtina=txtin;
   txtina{4}=wgy;
   txtina{1}=['Sign of Response -- ' txt1 '--' txt2 ' Model -- ' sitecode];
	txtina{1}=[txtina{1} ' (' txtyr ')'];
   datina{2}=flipud(G1);
   datina{4}=[1 64];
   figure(nfig2);
   colmap02(txtina,datina);
   if strcmp(kpost,'Yes'); % postscript file
      eval(['print -dpsc ' 'wind' int2str(nfig2) ';']);
   end
   
   
   %-----------  OE (1,1)
   
   Y1=Y(:,:,3);
   G1=G(:,:,3);
	txt2='OE(1,1)';
   nfig1=11; nfig2=12;
   
   txtin{1}=['Signal Strength -- ' txt1 '--' txt2 ' Model -- ' sitecode];
	txtin{1}=[txtin{1} ' (' txtyr ')'];
   datin{2}=flipud(Y1);
   figure(nfig1);
   colmap01(txtin,datin);
   if strcmp(kpost,'Yes'); % postscript file
      eval(['print -dpsc ' 'wind' int2str(nfig1) ';']);
   end

   datina = datin;
   txtina=txtin;
   txtina{4}=wgy;
   txtina{1}=['Sign of Response -- ' txt1 '--' txt2 ' Model -- ' sitecode];
	txtina{1}=[txtina{1} ' (' txtyr ')'];
   datina{2}=flipud(G1);
   datina{4}=[1 64];
   figure(nfig2);
   colmap02(txtina,datina);
   if strcmp(kpost,'Yes'); % postscript file
      eval(['print -dpsc ' 'wind' int2str(nfig2) ';']);
   end
end

disp('ALL DONE!');
   
   
  
	 

%******************  SUBFUNCTION TO DO THE OE MODELING
function [R,G]=respoe(C,D,endmo,kclim,kmode);

% A dedicated subfunction of resp02.m that either computes a
% a variance explained statistic for OE models with tree-ring 
% output and climate input, or computes the lag 1 acf
% of the climate input.  OE models (1,0), (2,0) and (1,1) are fit

% C is a col vector of the tree-ring series
% D is a lagged matric of climate series -- same years, as built in resp02
%		D has 37 columns
% endmo is ending month of the tree year (say, 8 for august)
% kclim is the the climate data type: 1==sum (e.g., ppt)
%		2== average (e.g., tmp)
% kmode is the mode of the run -- 1 is OE modeling,
%		2 is just get the lag1-1 autocorrelations of the climate series

nobs = size(C,1); % number of years for modeling

% Subtract mean from tree-ring series
C = C -mean(C);
cvar = var(C);  % variance of the output series


%----------  PREALLOCATE *********************************************

a=NaN;
R=a(ones(12,1),ones(12,1),ones(3,1)); % to hold variance explained
LG=R; % sign of largest irw + or -
Wzero=0;
W = Wzero(ones(12,1),ones(12,1),ones(3,1)); % will hold indicator of whether any of B or F params more than 
% 3 std devs from zero
g30=30;
G=g30(ones(12,1),ones(12,1),ones(3,1)); % indicator of whether one of following
% 15 -- explained variance >=10% and highest ir weight  -
% 30 -- explained variance less than 10%
% 45 -- explainde variance >=10% and higest ir weight +

J=[1:13:144]';

S5=zeros(nobs,12);

% Some quantities needed later for computing theoretical irf
ufake=zeros(nobs,1);
ufake(1)=1;


lastmo=37- (12-endmo);  % tie ending month to column of D

for k1=1:12;  % loop for each number of months in seasonal grouping
	top=lastmo-(k1-1):lastmo;  % rv.  Say [32 33] for k1 = 2,endmo=8
	S1 = top(ones(12,1),:); % dupe rows of top
	S2 = (0:11)';
	S3 = S2 (:,ones(k1,1)); % dupe cols of S2
	S4 = S1-S3;

	for k2 = 1:12;  % sum or ave the clim variable over months
		if k1==1; % special case for vector
			S5(:,k2) = D(:,S4(k2,:));
		else 
			S5(:,k2) = (sum((D(:,S4(k2,:)))'))';
			if kclim==2; % if a 'temperature like' variable
				S5(:,k2)=S5(:,k2)/k1; % convert from seasonal sum to seasonal ave
			else
				% no action needed, this is a 'precip like' variable
			end
		end ;  % of if loop
	end;  % of for loop k2

	% S5 is matrix of climate variable summed over months or averaged over months
	% S5 is specific to one value of number of months in season
	% cols of S5 are series with different ending months
	

	% Subtract means of seasonal climate series; then will have mean-zero series
	% for both climate and tree rings;  need this for OE modeling	
	S5mean = mean(S5);
	S5 = S5 - repmat(S5mean,nobs,1);


	% Decide if OE modeling or just compute first-order autocorrelation
	if kmode==1; % OE modeling
		% Loop over the climate series
		for n = 1:12;
         z=[C S5(:,n)];  % output in col 1, input in col 2

			 % I do not use this section now -- kept for reference
         % Compute impulse response function 
         %mlags=3;
         %arord=5;
         %[g,lag,sig99]=irf01(z,mlags,arord);
         %if g>=0; % largest of first mlags irf weights is positive
         %   LG(n,k1)=1;
         %else; % if largest weight negative (inverse relationship)
         %   LG(n,k1)=0;
         %end
        
       
			% Fit the OE models and associated info
			th1 = oe(z,[1 0 0]); 
			R(n,k1,1) = 1.0 - (th1(1,1)/cvar); %  1.0-(error variance/original variance
			yfake=idsim(ufake,th1);  % theoretical irf
			if yfake(1)>=0;
				g=1; % lag-0 irf weight positive
			else
				g=0;
			end
         LG(n,k1,1)=g;
         % Are any of the B or F parameters more than 3 std devs from zero?
         nparams=th1(1,5)+th1(1,8);
         if nparams~=1;  
            error('theta row 1 should give 1 param for OE(1,0) model');
         end
         % Get variances of estimate parameters; convert to standard
         % deviations; test for greater than one standard deviations
         % from zero
         cp = th1(3,1:nparams); % parameter estimates, B-operator, then F
         H = th1(4:(3+nparams),1:nparams); % variances of parameters
         c2 =1*sqrt(diag(H)); % one standard deviations of parameters
         wtsig = all(abs(cp)>abs(c2'));  % 1 if any of B or F params significant
         wtsig=logical(wtsig);
         if wtsig;
            W(n,k1,1)=1;
         end
      
         %---------OE(2,0) model
         th2 = oe(z,[2 0 0]);
			R(n,k1,2) = 1.0 - (th2(1,1)/cvar);
			yfake=idsim(ufake,th2);  % theoretical irf
			if (yfake(1)+yfake(2))>=0;
				g=1; % lag-0 irf weight positive
			else
				g=0;
			end
         LG(n,k1,2)=g;
         % Are any of the B or F parameters more than 3 std devs from zero?
         nparams=th2(1,5)+th2(1,8);
         if nparams~=2;  
            error('theta mtx row 1 should say 2 parmams for OE(2,0) model');
         end
         % Get variances of estimate parameters; convert to standard
         % deviations; test for greater than two standard deviations
         % from zero
         cp = th2(3,1:nparams); % parameter estimates, B-operator, then F
         H = th2(4:(3+nparams),1:nparams); % variances of parameters
         c2 =1*sqrt(diag(H)); % two standard deviations of parameters
         wtsig = all(abs(cp)>abs(c2'));  % 1 if any B and F params significant
         wtsig=logical(wtsig);
         if wtsig;
            W(n,k1,2)=1;
         end
         
         
         %---------  OE(1,1) model
         th3 = oe(z,[1 1 0]);
         R(n,k1,3) = 1.0 - (th3(1,1)/cvar);
			yfake=idsim(ufake,th3);  % theoretical irf
			% Find which impulse weight of first 5 is largest absolute value  
			[ybig,ibig]=sort(abs(yfake(1:5)));
			ythresh = 0.1 * abs(yfake(ibig(5)));
			ilarge=find(abs(yfake(1:5))>=ythresh);
			ilast = max([ibig(5) max(ilarge)]);
			gsum = sum(yfake(1:ilast));
			if gsum>=0;
				g=1;
			else
				g=0;
			end
         LG(n,k1,3)=g;
         % Are any of the B or F parameters more than 3 std devs from zero?
         nparams=th3(1,5)+th3(1,8);
         if nparams~=2;  
            error('theta mtx row 1 info should say 2 params for OE(1,1) model');
         end
         % Get variances of estimate parameters; convert to standard
         % deviations; test for greater than two standard deviations
         % from zero
         cp = th3(3,1:nparams); % parameter estimates, B-operator, then F
         H = th3(4:(3+nparams),1:nparams); % variances of parameters
         c2 =1*sqrt(diag(H)); % two standard deviations of parameters
         wtsig = all(abs(cp)>abs(c2'));  % 1 if any of B or F params significant
         wtsig=logical(wtsig);
         if n==6 & k1==1;
               disp('ditto'); % ditto
         end
        
         if wtsig;
            W(n,k1,3)=1;
         end
         
      end
      
       
	else; % kmode must be 2, and use only wants first-order autocorrel of clim
		S6=covf(S5,2);
		S6=S6(J,:);
		R(:,k1,1)= S6(:,2) ./ S6(:,1);
	
	end
   
   
end; % of loop for k1

%--------------- Fill G

% This obsolete version based 'signifificance' on explained variance 10% or mor
%for k5=1:3;
%	 LGsub = LG(:,:,k5);
%   Gsub=G(:,:,k5);
%   L15 = R(:,:,k5)>=0.10 & LGsub==0;
%   L45 = R(:,:,k5)>=0.10 & LGsub==1;
%   L15=logical(L15);
%   L45=logical(L45);
%   if any(any(L15));
%      Gsub(L15)=15;
%   end
%   if any(any(L45));
%      Gsub(L45)=45;
%   end
%   G(:,:,k5)=Gsub;
% end

% This version bases 'significance' on any of the estimated B or F coefs at
% least 1 std devs from zero and explained variance 5% or more
W = logical(W);
for k5=1:3;
   Rsub = R(:,:,k5);
	LGsub = LG(:,:,k5); % entries are 1 if + relation, 0 if negative
   Gsub=G(:,:,k5); % had been initialized as 30
   LW = W(:,:,k5); % logical, 1 if any signif B or F weights
   L15 = LW & Rsub>=0.25 & LGsub==0; % relation strong and direction of relation negative
   L45 = LW & Rsub>=0.25 & LGsub==1; % relation strong  ... postive
   % Recall that G initialized as 30.  Replace selected entried of G as needed
   if any(any(L15));
      Gsub(L15)=15;
   end
   if any(any(L45));
      Gsub(L45)=45;
   end
   G(:,:,k5)=Gsub;
 end

  


