function datout=pdsi1(datmon,datother,snowinf,textin,kopt,penopts,datpen)
% pdsi1:  Palmer Drought Severity Index
% CALL: datout=pdsi1(datmon,datother,snowinf,textin,kopt,penopts,datpen);
%
%*******************  IN **************************
%
% datmon {} monthly 13-col input data series
%	{1} P (? x 13)r precip (in)
%	{2} T (? x 13)r temperature (oF)
%
% datother {} additional required input
%	{1} lat (1 x 1)r latitude (decimal degrees)
%	{2} watcap(1 x 2)r soil moisture info
%		awcs (1 x 1)r available water capacity, surface layer (in)
%		awcu (1 x 1)r available water capac., underlying layer (in)
%	{3} YRS (2 x 2)i  start and end years
%		row 1: desired pdsi output
%		row 2: period for normals & "CAFEC" computations
%	{4} dayz lookup table of pctg of possible sunshine
%
% snowinf{} input relating to snow storage model -- see notes
%		Tcrit (1 x 1)r critical monthly mean temperature (oF) below which
%			ppt is assumed to be snow.  Typically 30.2 oF.  May set to []
%			if kopt(1)==1.
%		mons1 (1 x 1)i start, end months of 'snow period', typically [10 5], 
%			meaning Oct-May.  This is the group of months snow is allowed to
%			fall in.  Then any snow that accumulated must be totally melted
%			by the start of the next snow year.  For the [10 5] example,
%			snow must be melted by the end of september. User flagged with
%			error if snow falls later than mons(2) month. May set to []
%			if kopt(1)==1;
%		melttbl (? x ?)r  lookup table to be developed relating number of
%			months for the melt to total snow accumulation in cold preceding
%			months and T of the first month with T>=Tcrit. May set to [] if
%			kopt(1)==1;
% textin {} text information
%	{1} (1 x ?)s stnname: name of station, region, or whatever -- used in plot titles
%
% kopt (1 x ?)i   general options
%		kopt(1) method for PE computation
%			==1 Thornthwaite
%			==2 Penman
%			==3 Scaled pan evaporation
%		kopt(2) mode of run
%			==1 routine (sparse interaction and graphics)
%			==2 instructional or exploratory mode
%		kopt(3) snow option
%			==1 treat all ppt as rain, regardless of temperature
%			==2 redistribute monthly ppt if T<Tcrit
% penopts (1 x ?)i   options on Penman input data
%		penopts(1) wind
%			==1 monthly time series input
%			==2 12 long-term monthly means only
%		penopts(2) relative humidity
%			==1 monthly time series input
%			==2 12 long-term monthly means only
%
% datpen {} Penman input data  ;  see notes
%	{1} wind speed
%	{2} relative humidity
%
%
%
%************************ OUT *************************
%
% datout {} output data
%  {1} Z (? x 13)r "Z" index
%  {2} X (? x 13)r Palmer Index
%  (3) XM (? x 13)r Modified Palmer Index
%  {4} W (? x 13)r Ave monthly soil moisture (in)
%  {5} RO (? x 13)r model runoff
%  {6} S1 (? x 13 x 10)r effective ppt  = max(0,p-f*pe), f=.1,.2,...1.0
%     The value for a given month/yr is the precip minus a fraction of the
%     estimated potential evapotranspiration.  If the difference is less than 
%     zero, the effective ppt for the  month is set to zero 
%  {7} P (? x 13)r pcp data, same size as Z, X, etc
%  {8} T {? x 13}r tmp data, same size ...
%**********************  NOTES ************************
%
% snowinf{} empty is OK if kopt(3)==1
% 
% Normals period assumed to be enclosed within output period
% In other words, YRS(2,:) enclosed within YRS(1,:)
%
%
% penopts, datpen can be set to [] if not using penman method
%
% Snow storage model optionally redistributes monthly ppt for those years with
% snow accumulation.  Snow is assumed if the monthly T <Tcrit. Accumulated
% snow is assumed to melt completely in the first  1-several months with
% T>=Tcrit.  The number of months (nmelt) is interpolated from a 2-d
% lookup table based on the total snow accumulation for the year and 
% the temperature for the month following the last month with 
% T<Tcrit.  So far, this is a constant lookup table, because I 
% have not developed an empirical relation ship. I assume a 2-month melting
% period.
%
% The percentage of accumulated snow assigned as redistributed ppt to the
% first warm months is proportional to T-Tcrit for those months.  For
% example, if the first two warm months have T = 32.2 and 40.2, the
% weights are 2/10 and 8/10
%
%--------------------
% Steps 
%
%	Check input
% Make pointers to desired years of calib period and output
% Compute PE according to selected PE option 
%		(Heat Index for Thornthwaite PE should be based on calib-pd normals)
% Redistribute montly P, if needed, according to selected snow option
% Get calib-period PE, P, T and compute climatic 'normals'
% Build synthetic time series of the normals of PE and P, repeated 10 yr
% Call soilmoi1.m using saturated initial conditions and the synthetic
%		10-yr seriesto get revised initial Jan starting soil moistures
% Call soilmoi1.m using the actual calib-period P,PE to get quantities
%		for  coefficients  of evapotranspiration, recharge, runoff, loss 
% Compute the normals and the coefficients
% Call soilmoi1.m on full-length actual data to get suite of accounting
%		time series: PR, PRO, etc
% Apply the coefficients from earlier step to these time series to 
%		get CAFEC precipitation and excesses and deficiencies (eqn 15, p. 15)
% Compute the 12 monthly weighting factors K' (eqn 26, p. 25)
%		(Be sure to base this on the average D, etc for calib period, not the
%		full period)
% Compute the 12 monthly values for K (eqn 27, p. 26)
% Compute z values from K and d
% Compute PDSI values X
% Gather and store key outputs into cell variable
% Return
%
%*************************  FUNCTIONS CALLED
%
% pdsix.m -- PDSI "X" values, from  Z index
% pethorn.m -- compute PE from monthly mean T
% snowppt.m <optional>utility function to redistribute monthly ppt over months by
%		adjusting for storage of ppt in a "snow buffer". Required inputs are
%		same-size matrices T and P.  snowppt.m returns the revised P.  Function
%		snowppt.m is not yet written.  Is needed only if a monthly T is below
%		Tcrit, and if kopt(3)==2
%
%------ highrow.m -- utility function to find highest row with nonzero 
%			logical element in each col of a matrix 
%
% soilmoi1.m -- soil moisture accounting
% 


%----------------  Prelims
close all
a=NaN;



% ----------------- CHECK INPUT DIMENSIONS

% datmon
[mtemp,ntemp]=size(datmon);
if mtemp~=1 & ntemp~=2,
	error('datmon should be row-cell of length 2');
end

% datother
[mtemp,ntemp]=size(datother);
if mtemp~=1 & ntemp~=4,
	error('datother should be row-cell of length 4');
end

% textin
[mtemp,ntemp]=size(textin);
if mtemp~=1 & ntemp~=1,
	error('datother should be row-cell of length 1');
end

% kopt
[mtemp,ntemp]=size(kopt);
if mtemp~=1 & ntemp~=3,
	error('kopt should be rv of length 3');
end

% snowinf
if kopt(3)==1;
	snowinf=[];
else
	[mtemp,ntemp]=size(snowinf);
	if mtemp~=1 | ntemp~=3;
		error('snowinf should be row-cell of length 3');
	Tcrit=snowinf{1};
   [m1temp,n1temp]=size(Tcrit);
   end
	if ~ (m1temp==1 & n1temp==1);
		error('Tcrit should be scalar');
	end
	mons1=snowinf{2};
	[m1temp,n1temp]=size(mons1);
	if ~(mons1(1)==1 & mons1(2)==2);
		error('mons1 must be 1 x 2');
	end
	melttbl=snowinf{3};
	error('hey, I have not built this yet');
end	


% penopts and datpen
if kopt(1)==1; % using thornthwaite method for pe
	penopts=[];
	datpen=[];
else
	[mtemp,ntemp]=size(penopts);
	if mtemp~=1 & ntemp~=2,
		error('penopts should be rv of length 2');
	end
	% datpen
	[mtemp,ntemp]=size(datpen);
	if mtemp~=1 | ntemp~=2,
		error('datpen should be row-cell of length 2');
	end
end


%--------------- Unpack YRS 

YRS = datother{3};
if size(YRS,1)~=2 | size(YRS,2)~=2,
	error('YRS must be 2 x 2');
end
	

%----------------  Check YRS against P, T

P=datmon{1}; % monthly precip data
T=datmon{2}; % monthly tmp data
[mtemp,ntemp]=size(P);
if ntemp~=13;
	error('P must be col size 13');
end
[mtemp,ntemp]=size(T);
if ntemp~=13;
	error('P must be col size 13');
end

%-----------------------------------
% Find the longest common period of P and T.  This will define the
% start and end of the hydrologic accounting period

% col vectors of years for P, T
yrpa=P(:,1);
yrta=T(:,1);

% Start and end years of common period for P, T
yrgo1=max(min(yrpa),min(yrta));
yrsp1=min(max(yrpa),max(yrta));


yr1=(yrgo1:yrsp1)';  % year vector for common perod
if YRS(1,1)<yrgo1 | YRS(1,2)>yrsp1;
	error('Desired output years YRS(1,:) inconsistent with common period for P,T');
end
if YRS(2,1)<yrgo1 | YRS(2,2)>yrsp1;
	error('Desired CAFEC years YRS(2,:) inconsistent with  common period for P,T');
end	

%---------------   station name, latitude and water capacities
lat=datother{1}; % latitude of station or region, assumed dec degrees
if lat<0 | lat>90;
	error('lat must be between 0 and 90');
end
watcap=datother{2}; % avail water capacities
[mtemp,ntemp]=size(watcap);
if ~(mtemp==1 & ntemp==2);
	error('watcap must be row-cell of length 2');
end
awcs=watcap(1);  % avail wat capac, sfc
awcu=watcap(2);  % avail wat capac, underlying

%-------------   Unpack lookup table for pctg possible sunshine
dayz=datother{4};
[mtemp,ntemp]=size(dayz);
if ~(mtemp==51  & ntemp==12);
	error('lookup table dayz must be 51 x 12');
end


%-----------------  Unpack text info
stnname=textin{1}; % station name
if ~isstr(stnname);
	error('stnname should be string');
end

clc;
disp(['STARTING ANALYSIS OF: ' stnname]);


%**************** MAKE POINTERS TO YEARS IN P, T ****************

% To common years for input P,T.  This is the 'full' period for computation of PDSI
% Logical row pointer is relative to input years of P, T
LP1 = yrpa>=yrgo1 & yrpa<=yrsp1;  % to common years  for P
LT1 = yrta>=yrgo1 & yrta<=yrsp1;  % to common years  for t

% Cull full-period P, T;  P, T will now cover same years
P=P(LP1,:);
T=T(LT1,:);
nyrsfull=size(P,1);


% Pointer to full-period P, T indicating years for 
% calibration period, also used for CAFEC coefficients and
% normals
LPT2 = yr1>=YRS(2,1) & yr1<=YRS(2,2);
yr2=yr1(LPT2);
nyrs2=length(yr2);


% To desired output years 
LPT1 = yr1>=YRS(1,1) & yr1<=YRS(1,2);
yrout=yr1(LPT1);


%************************ COMPUTE PE ****************************

disp('   Compute PE');

if kopt(1)==1; % Thornthwaite method 
	% Compute PE, using calib period for Heat Index
	% N means N hemis
	PE=pethorn(T,YRS(2,:),lat,dayz,'N');
else
	error('Only Thornth PE method programmed so far');
end

% PE is now the 13-col matrix of full-period PE



%****************** REDISTRIBUTE P ACCORDING TO SNOW OPTION

if kopt(3)==1;
	% No action necessary. Treat snow as rain.
	Porig=P; % Porig is same as P in this branch
else;
   disp('   Snow option redistribution of monthly precip');
	error('Snow option redistrib of monthly P not yet programmed');
	Porig=P; % store precip before redistribution
	% This option would give revised monthly P
end


%******  GET CALIBRATION PERIOD P,T,PE AND COMPUTE NORMALS

disp('   Soil Moisture Accounting -- pass 1 to get starting Jan soil moistures');

Pbar = mean(P(LPT2,2:13));
Tbar = mean(T(LPT2,2:13));
PEbar= mean(PE(LPT2,2:13));


%*************** BUILD SYNTHETIC 10-YR TIME SERIES OF NORMALS OF P, PE

Psyn=[(1:10)' repmat(Pbar,10,1)];
PEsyn = [(1:10)'  repmat(PEbar,10,1)];


%**************  CALL SOILMOI1.M TO GET INITIAL SOIL MOISTURE FOR LATER RUN
%
% Strategy here is to run the soil moisture balance for 10 years, starting
% with conditions saturated to start the first january.  By the last, or 10th
% January, will have arrived at more reasonable values for jan 1 soil moisture

% String the 12-col data matrices to col vectors of precip and potential evap
ptemp= (Psyn(:,2:13))';
ptemp=ptemp(:);

petemp=(PEsyn(:,2:13))';
petemp=petemp(:);

% Soil moisture accounting, assuming saturated conditions jan 1 of first year
meat=soilmoi1(ptemp,petemp,awcs,awcu,awcs,awcu);


% Note that meat{4} has start of month soil moisture in sfc layer
%  and that meat{5} has it for underlying layer
% For a 10-year (48  month) run, the starting jan soil moistures for the
% last year are in element 37 of 48.
ssgo=meat{4};
sugo=meat{5};
SSGO=(reshape(ssgo,12,10))';
SUGO=(reshape(sugo,12,10))';

ssgo=SSGO(10,1); % starting jan soil moisture in sfc layer
sugo=SUGO(10,1); % starting jan soil moist in underlying layer
clear meat;


%*********  SOIL MOISTURE ACCOUNTING USING CALIB PERIOD TIME SERIES
% DATA AND PREVIOUSLY ESTIMATED STARTING JAN SOIL MOISTURES

disp('   Soil Moisture Accounting -- Pass 2 to get quantities for CAFEC coefs');

% Get matrices of calib-period P, PE, and convert to col vect
ptemp=(P(LPT2,2:13))';
ptemp=ptemp(:);
petemp=(PE(LPT2,2:13))';
petemp=petemp(:);

% Soil moisture accounting
meat=soilmoi1(ptemp,petemp,awcs,awcu,ssgo,sugo);


%**************  COMPUTE NORMALS AND CAFEC COEFS

% recall that nyrs2 is number of years in calib period

et=meat{17};
ET=(reshape(et,12,nyrs2))';

% Recharge and potential recharge
r=meat{11};
R = (reshape(r,12,nyrs2))';
pr=meat{12};
PR = (reshape(pr,12,nyrs2))';

% Runoff and potential runoff
ro=meat{13};
RO = (reshape(ro,12,nyrs2))';
pro=meat{14};
PRO = (reshape(pro,12,nyrs2))';

% Loss and potential loss
loss=meat{15};
LOSS = (reshape(loss,12,nyrs2))';
ploss=meat{16};
PLOSS = (reshape(ploss,12,nyrs2))';

ETbar=mean(ET);
Rbar=mean(R);
PRbar=mean(PR);
RObar=mean(RO);
LOSSbar=mean(LOSS);
PRObar=mean(PRO);
PLOSSbar=mean(PLOSS);


% alpha
anan=NaN;
alpha=anan(:, ones(12,1)); % initialize as NaN
L1=PEbar==0 & ETbar==0;
L2=PEbar==0 & ETbar~=0;
L3=~ (L1 | L2);
if any(L1);
	alpha(L1) = 1;
end
if any(L2);
	alpha(L2) = 0;
end
if any(L3);
	alpha(L3)=ETbar(L3) ./ PEbar(L3);
end

% beta
beta=anan(:, ones(12,1)); % initialize as NaN
L1=PRbar==0 & Rbar==0;
L2=PRbar==0 & Rbar~=0;
L3=~ (L1 | L2);
if any(L1);
	beta(L1) = 1;
end
if any(L2);
	beta(L2) = 0;
end
if any(L3);
	beta(L3)=Rbar(L3) ./ PRbar(L3);
end


% gamma
gamma=anan(:, ones(12,1)); % initialize as NaN
L1=PRObar==0 & RObar==0;
L2=PRObar==0 & RObar~=0;
L3=~ (L1 | L2);
if any(L1);
	gamma(L1) = 1;
end
if any(L2);
	gamma(L2) = 0;
end
if any(L3);
	gamma(L3)=RObar(L3) ./ PRObar(L3);
end


% delta
delta=anan(:, ones(12,1)); % initialize as NaN
L1=PLOSSbar==0;
if any(L1);
	delta(L1)=0;
end
L2=~L1;
if any(L2);
	delta(L2) =   LOSSbar(L2) ./ PLOSSbar(L2);
end


%********* SOIL MOISTURE ACCOUNTING ON FULL-LENGTH TIME SERIES

disp('   Soil Moisture Accounting -- Pass 3 on full length time series');

% Get matrices of full period P, PE, and convert to col vect
ptemp=(P(:,2:13))';
ptemp=ptemp(:);
petemp=(PE(:,2:13))';
petemp=petemp(:);

% Soil moisture accounting
meat=soilmoi1(ptemp,petemp,awcs,awcu,ssgo,sugo);

et=meat{17};
ET=(reshape(et,12,nyrsfull))';

% Recharge and potential recharge
r=meat{11};
R = (reshape(r,12,nyrsfull))';
pr=meat{12};
PR = (reshape(pr,12,nyrsfull))';

% Runoff and potential runoff
ro=meat{13};
RO = (reshape(ro,12,nyrsfull))';
pro=meat{14};
PRO = (reshape(pro,12,nyrsfull))';

% Loss and potential loss
loss=meat{15};
LOSS = (reshape(loss,12,nyrsfull))';
ploss=meat{16};
PLOSS = (reshape(ploss,12,nyrsfull))';


% Ave of start of month and end of month soil moisture, combined sfc and underlying
smave=meat{10};
W=(reshape(smave,12,nyrsfull))';
W=[yr1 W];

% Strip year column off PE, P, T, to make compatible for matrix operations
% Note that to restore, just need to slap on yr1 in as the first col
PE(:,1)=[];
P(:,1)=[];
T(:,1)=[];

% Effective precip, defined as monthly precip minus a fraction of the PE
S1 = a(ones(nyrsfull,1),ones(13,1),ones(10,1));
for ncount = 1:10;
   nfract = ncount/10;
   S1(:,:,ncount) = [yr1 max(P- nfract * PE,0)];
end

   


%************  COMPUTE CAFEC PRECIP -- see p. 14, eqns 10-14

EThat=(repmat(alpha,nyrsfull,1)) .* PE;
Rhat = (repmat(beta,nyrsfull,1)) .* PR;
ROhat = (repmat(gamma,nyrsfull,1)) .* PRO;
LOSShat=(repmat(delta,nyrsfull,1)) .* PLOSS;
Phat = EThat + Rhat + ROhat - LOSShat;



%**************** COMPUTE PRECIP EXCESSES AND DEFICIENCIES
disp('   Compute precip departures and K factors');

d = P -Phat;


%********************* THE 12 MONTHLY K' FACTORS ( P. 25, EQN 26)

% Get the departures for the calibration period
dcalib = d(LPT2,:);

% Get the mean absolute departures for each month of year
Dbar= mean(abs(dcalib));

% Compute K' values from eqn 26
Kprime =     1.5*log10((((PEbar+Rbar+RObar)/(Pbar+LOSSbar)) + 2.80) ./ Dbar) +0.50;



%*********************** FINAL ADJUSTMENT OF THE MONTHLY K' VALUES

% Compute weighted ave departures for each month
DKprime= Dbar .* Kprime;

% Compute the sum of these 12 values
denom = sum(DKprime); % denominator of eqn 27, p. 26


% Compute final 12 K values
K = (17.67/denom) * Kprime;



%******************  COMPUTE Z INDICES FROM d and K
disp('   Compute Z index');

% d is a time series matrix (something x 12) and K is 12 monthly constants
% Need to dupe K to size of d before multiplying
Z = (repmat(K,nyrsfull,1)) .* d;


%***************  COMPUTE PDSI FROM THE Z INDICES
disp('   Compute PDSI')
% String z into a col vector first
Ztemp=Z';
z=Ztemp(:);
[x,xm] = pdsix(z);
X = (reshape(x,12,nyrsfull))';
XM = (reshape(xm,12,nyrsfull))';

datout={[yr1 X],[yr1 XM],[yr1 Z], W,[yr1 RO], S1,[yr1 P],[yr1,T]};

