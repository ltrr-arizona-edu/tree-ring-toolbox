% screen3.m -- a script file
%
% Third in a sequence of functions for screening chrons vs
% climatic series. 
%  screen1.m -- chrons vs nearest ns stations 
%  screen2.m -- uses results from screen1.m to make pointers to 
%		"best" station to pair with each site.  "Best" by length
%		of overlap data for modeling and proximity
%  screen3.m -- model chrons vs their paired climate series as
%		identified in screen2.m
%  screen4.m -- prepare results for mapping by surfer
%
% Originally written as a function, but because of "stack" errors,
% used as a script file. The input "arguments" listed below must
% be in the base workspace.  The output "arguments" can be culled
% and saved by the command in the notes section below.
%
% D Meko 6-28-96
%
%************* INPUT - from screen-prompted input or files
%
% T (mT x nT)r chronology indices; NaN filled; variable start and
%		end years allowed; mT years, nT chronologies
% C (mC x nC)r climate variable; NaN filled; variable start and 
%		end years allowed; mC years, nC stations
% YRS (2 x 2)i start, end years of 
%		row-1: T
%		row-2: C
% mx (1 x 5)i maximum setting for 
%  1-number of stations to find in search radius
%  2-number of B-operator parameters in a model
%  3-number of F-operator parameters
%  4-number of total parameters in model
%  5-number of ranked models to consider (nranked)
% O (1 x 1)i options -- userprompted
%		1- 1=regular mode;  2=diagnostics mode
%
%---- next info normally loaded from distinfa.mat ---
%
% N (nT x 1)i number of stations found in search radius;
%		nT chronologies
% W (nT x mxs)i which stations (columns of C) grouped with this
%		chronology
% D (mD x nD)r distance (km) from chronology to each station in
%		search radius;  mD==nT; nD==mxs
%
%---- next info normally loaded from screen2.m output file,
%	named something like AP1OS2.MAT
%
% w (mw x 1)i stations to use
% iw (mi x 1)i first(1) or second(2) class pair
% 
%************************* OUT ARGS **********************
%
%
% K2 (nT x 1)i  significant signal?  meaning are any of the
%		parameters of the full-period model signifantly different
% 		from zero (2+ sdevs); 1=yes 0=no
% OB (nT x 1)i  order of B-operator in model
% OF (nT x 1)i  order of F-operator in model
% R1 (nT x 1)r  corr coefficient between y(t) and u(t)
% R2 (nT x 1)r  corr coefficient between y(t) and [B/F] u(t)
% G (nT x mxs)i number of significant (99% level) autocorrelation
%		coefficients  of residuals of model
% H (nT x 1)i number of significant (99% level) cross-
%		correlations between input u(t) and residuals e(t)
% S (nT x 1)r variance-explained fraction for best OE model
%		Computed as 1 - (var(residuals)/var(y)) 
%		for final model fit to entire data period
% NY (nT x 1)i number of years used to fit final OE model
%
%****************** STEPS IN SCREENING*************************
%
%  For tree-set "a", "p"recipt, season "1", data files are
%
% Check input data
% load ap1is1.mat
%    Change O if necessary
% load distinfa.mat
% load c1104a.mat 
% load ap1wa1.mat
% Model each chronology -- output series y(t)
%	Pull out submatrices of columns
%	Pull row subset for valid overlap years
%		Loop over potential models
%			Find best model
%		End Loop
%		Fill slots in R1, R2, etc
%	End Loop
% 
%
%******************** NOTES *************************************
%
% Function oesplit1.m is called to perform split-sample modeling;
% calibrate on first half and verify on second; then vice versa;
%
% Details of the split-sample modeling are covered in 
% oesplit1.m
%
% Modeling period is flexible to differing time coverage of the
% climate station series and chronologies;  model is fit to 
% full common period for any pair
%
% The call to near3 is hard-coded to get the nearest mx(1) stations,
% regardless of whether they are in the search radius.  To change
% this, see instructions for near3.m and pass a different option
% to near3.m
%
% Notice the loop n=1:nT.  Normally you run this function on all 
% chronologies.  But you can control it by specifying some
% continuous sub-group of chronologies,by hard coding a change, like
%
% for n = 45:45 ... to  just model chron 45
%
% A combination of this control in the n loop and setting O(1)==2
% can be used to get detailed graphics model results in the figure
% windows for a specific chron/station pair.
%
%****************************************************

% Help in file selection by specifying letter code for tree-ring set
char=input('Letter code for tree-ring matrix: ','s');


% Get the tree-ring data tsm
fn5=['tsm' char '.mat'];
txt5='TSM? (in D:\wrk6)--  tree-ring indices';
[file5,path5]=uigetfile(fn5,txt5);
eval(['load ',path5,file5]);  % tree-data will be in X, year in col1
X(:,1)=[]; % get rid of year
T=X; 
clear X N


% Get the climate data tsm
fn5=['C*.mat'];
txt5='CP#.mat  -- climate data,season #';
[file5,path5]=uigetfile(fn5,txt5);
eval(['load ',path5,file5]);  % climate data is C

% T and C should now be in workspace

% Get the years info and miscellanous settings
fn5=[char '*.mat'];
txt5=[char 'P#IS1.MAT -- misc screen1.m input']; 
[file5,path5]=uigetfile(fn5,txt5);
eval(['load ',path5,file5]); 

% Get the site-to-station distance info
fn5=['DIST' char '.mat'];
txt5='DIST?.MAT -- distance info'; 
[file5,path5]=uigetfile(fn5,txt5);
eval(['load ',path5,file5]); 

% Get the "class-of-model" info as produced by screen2.m
fn5=[char '*.mat'];
txt5=[char 'P#OS2.MAT -- screen2.m output']; 
[file5,path5]=uigetfile(fn5,txt5);
eval(['load ',path5,file5]); 


clear path5 file5
pack

% User prompt for whether quick or diagnostics mode
O=1;
O=input('Quick mode (1) or diagnostics mode (2): ');


a=NaN;


%****************** CHECK INPUT DATA, SIZE, ALLOCATE

[mT,nT]=size(T);
[mC,nC]=size(C);

[m1,n1]=size(O);
if m1~=1 | n1~=1
	error('O must be 1 x 1')
end
if O(1) ~=1 & O(1)~=2,
	error('Allowable O(1) values are 1 and 2')
end

[m1,n1]=size(YRS);
if m1~=2 | n1~=2
	error('YRS must be 2 x 2');	
end


[m1,n1]=size(mx);
if m1~=1 | n1~=5,
	error('mx must be 1 x 5')
end
maxs=mx(1); % maximum number of stations to grab within search
	% radius; use the nearest maxs stations
nb = mx(2); % maximum allowable order of B operator
nf = mx(3); % maximum allowable order of F operator
ntot = mx(4); % maximum total number of parameters allowed
nranked=mx(5); % number of "good" models to consider before
		% selecting best OE model


% Initialize some matrices
R1=a(ones(nT,1),ones(maxs,1)); % Correlations y(t) vs u(t)
R2=a(ones(nT,1),ones(maxs,1)); % Correlations y(t) vs Bu(t)
S=a(ones(nT,1),ones(maxs,1)); % Variance-explained fraction
Smore=a(ones(nT,1),ones(maxs,1)); 
OB=a(ones(nT,1),ones(maxs,1)); % B-operator order for best model
OF=a(ones(nT,1),ones(maxs,1)); % F-operator order for best model
K2=a(ones(nT,1),ones(maxs,1)); % any signif params in final model?
H= a(ones(nT,1),ones(maxs,1)); % number of sig autoc of residual
G= a(ones(nT,1),ones(maxs,1)); % number of sig cc of res with input
Q=a(ones(nT,1),ones(maxs,1)); % "way" best model was finally selected
NY=a(ones(nT,1),ones(maxs,1)); % number of years for final OE model


%**************** ADJUST CLIMATE MATRIX C SO ONLY
% THE LONGEST CONSEC NON-NAN SEQUENCES REMAIN NON-NAN
% In other words, the result will have no internal imbedded NaNs
for n = 1:nC;
	c1=C(:,n); % temporary storage of this stations climate vector
	c2=a(ones(mC,1),:); % same size vector of NaN
	[ii1,ii2]=consec(c1); % start,end row indices of good stretch
	nsub = ii2-ii1+1; % number of values in good stretch
	c2(ii1:ii2)=c1(ii1:ii2); % insert the good stretch in c2
	C(:,n)=c2; % replace column of C
end
clear c1 c2 nsub	



%**************** GET NEAR-STATIONS INFORMATION **************
%
% This step now down outside of screen3.m beforehand by steps listed
% below:
%
%k1=[1 2];
%PT= dec2dms2(TC,k1); % change coordinates from map long-lat to regular
%PC= dec2dms2(CC,k1);
%P=gcdist(PT,PC); % calc dist (km) from each tree-site to each station
%[D,W,N,I2]=near3(P,ds,maxs,maxs,1);
% D, W, N defined before;  I2 is pointer to valid chronologies, 
%		in other word, chronologies with at least one climate station
%		in the search radius
%
% Typically, would have saved
% D W N I2 in file, say, distinfa.mat
%
%*************** LOOP OVER CHRONOLOGIES ************************


pack


% Build a matrix of OE model structures.  These will be the
% candidate models.
nn=struc3(1:nb,0:nf,0,ntot); % no-delay models only

% Start the loop over chronologies

% If diagnostics mode, do only one chronology, otherwise
% will do all unless you manually change ngo and nsp
% to restrict it -- to say "ngo=25, nsp=30"
if O(1)==2; % diagnostics, one-chronology mode
	nthis = input('Chronology # to look at');
	ngo=nthis;
	nsp=nthis;
elseif O(1)==1;
	ngo=1;
	nsp=nT
else
	error('Invalid O(1) value')
end

for n =ngo:nsp; %1:nT;
	ns=N(n); % number of stations found for this chronology
	C1=C(:,W(n,1:ns)); % col sub-matrix of the "near" stations
	z= T(:,n); % all rows of data for this chronology
	c=C(:,w(n)); % all rows of the station climate series	
	disp(['On chron ',int2str(n)])

	% Get the valid data for common period; need input climate
	% series u and output chronology y, and year vector yrs
	[u,y,yrs]=overlap(c,z,flipud(YRS));

	% Get correlation coefficient between u and y
	rr1= corrcoef([y u]);
	r1=rr1(1,2);

	nyrs = length(u); % number of years for OE modeling

	%  OE-model the series with split-sample validation
	Z=[y u];
	%[s,r2,ob,of,f1,f2,way,k2]=oesplit1(Z,yrs,nn,nranked,O(1));
   % O(1)==1 means speedy mode; ==2 means diagnostic mode
   [s,r2,ob,of,f1,f2,way,k2,smore]=oefit(Z,yrs,nn,nranked,O(1));


	% Store the results for this chronology/station
	K2(n)=k2;
	OB(n)=ob;
	OF(n)=of;
	G(n)=f1;
	H(n)=f2;
	R1(n)=r1;
	R2(n)=r2;
   S(n)=s;
   Smore(n)=smore;
	NY(n)=nyrs;
   Q(n)=NaN;
   %Q(n)=way;
end; % of "n=" loop over chronologies1


set2=' K2 OB OF R1 R2 G H S NY ';

% Store result in a .mat file of the user's choice
disp('Example file name for next prompt is ap1os3.mat')
disp(' meaning tree-ring set ''a'', data precip')
disp(' season 1 -- where 1 is nov-apr')
disp(' o  -- output')
disp(' s3 -- from program screen3.m')
disp(' press return')
pause

fn1=[char '*.mat'];
txt1='?P#OS3.mat -- output from screen3.m';
[file1, path1] = uiputfile(fn1,txt1);
eval(['save ',path1,file1,set2])
