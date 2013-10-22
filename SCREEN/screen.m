% screen.m -- a script file
%
% Screen network tree-ring chronologies against climate;
% Originally written as a function, but because of "stack" errors,
% used as a script file. The input "arguments" listed below must
% be in the base workspace.  The output "arguments" can be culled
% and saved by the command in the notes section below.
%
% D Meko 4-30-96
%
%************************* IN ARGS ************************
%
% T (mT x nT)r chronology indices; NaN filled; variable start and
%		end years allowed; mT years, nT chronologies
% C (mC x nC)r climate variable; NaN filled; variable start and 
%		end years allowed; mC years, nC stations
% YRS (2 x 2)i start, end years of 
%		row-1: T
%		row-2: C
% TC (mTC x nTC)i decimal-degree long-lat file for T
% CC (mCC x nCC)i decimal-degree long-lat file for C
% ds (1 x 1)r search radius (km)
% mx (1 x 5)i maximum setting for 
%  1-number of stations to find in search radius
%  2-number of B-operator parameters in a model
%  3-number of F-operator parameters
%  4-number of total parameters in model
%  5-number of ranked models to consider (nranked)
% O (1 x 2)i options
%		1- 1=regular mode;  2=diagnostics mode
%		2- 1=model all nearby stations; 2=single-station modeling
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
%---- next info normally loaded from wa1.mat, or something like it
%
% w (mw x 1)i stations to use
% iw (mi x 1)i first(1) or second(2) class pair
% 
%************************* OUT ARGS **********************
%
%
% K2 (nT x mxs)i  significant signal?  meaning are any of the
%		parameters of the full-period model signifantly different
% 		from zero (2+ sdevs); 1=yes 0=no
% OB (nT x mxs)i  order of B-operator in model
% OF (nT x mxs)i  order of F-operator in model
% R1 (nT x mxs)r  corr coefficient between y(t) and u(t)
% R2 (nT x mxs)r  corr coefficient between y(t) and B u(t)
% G (nT x mxs)i number of significant (99% level) autocorrelation
%		coefficients  of residuals of model
% H (nT x mxs)i number of significant (99% level) cross-
%		correlations between input u(t) and residuals e(t)
% S (nT x mxs)r variance-explained fraction for best OE model
%		Computed as 1 - (var(residuals)/var(y)) 
%		for final model fit to entire data period
% NY (nT x mxs)i number of years used to fit final OE model
%
%****************** STEPS IN SCREENING*************************
%
% Check input data
% Get distance settings, sites in search radius
% Model each chronology -- output series y(t)
%	Pull out submatrices of columns
%	Loop over "near" climate stations -- input series u(t)
%		Pull row subset for valid overlap years
%		Loop over potential models
%			Find best model
%		End Loop
%		Fill slots in R1, R2, etc
%	End Loop
% 
%
%******************** NOTES *************************************
%
% Function oesplit.m is called to perform split-sample modeling;
% calibrate on first half and verify on second; then vice versa;
%
% Details of the split-sample modeling are covered in 
% oesplit.m
%
% Coordinate file TC and CC are assumed to be in mapping format,
% that is, decimal long (negative west) followed by decimal
% latitude
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
% To save output described above, do this
% save hey N W D OB OF R1 R2 G H S NY
%
% I first tried writing this code as a function, but the function
% bombed out with a nebulous stack error.  Maybe I had too many
% input or output arguments.  If so, I could possibly re-code by
% combining some matrices into larger matrices. Oh well, just in
% cas, here was the function call
%function [N,W,D,OB,OF,R1,R2,G,H,S,NY,Q]=screeny(T,C,YRS,TC,CC,ds,mx,O)
%[N,W,D,OB,OF,R1,R2,G,H,S,NY,Q]=screeny(T,C,YRS,TC,CC,ds,mx,O);
%
%****************** CHECK INPUT DATA, SIZE, ALLOCATE

[mT,nT]=size(T);
[mC,nC]=size(C);
[m1,n1]=size(TC);
if m1~=nT | n1~=2;
	error('TC must have 2 cols and be same row size as T');
end
[m1,n1]=size(CC);
if m1~=nC | n1~=2,
	error('CC must have 2 columns and be same row-size as C')
end
[m1,n1]=size(O);
if m1~=1 | n1~=2
	error('O must be 1 x 2')
end
if O(1) ~=1 & O(1)~=2,
	error('Allowable O(1) values are 1 and 2')
end
if O(2) ~=1 & O(2)~=2,
	error('Allowable O(2) values are 1 and 2')
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

a=NaN;

% Initialize some matrices
R1=a(ones(nT,1),ones(maxs,1)); % Correlations y(t) vs u(t)
R2=a(ones(nT,1),ones(maxs,1)); % Correlations y(t) vs Bu(t)
S=a(ones(nT,1),ones(maxs,1)); % Variance-explained fraction
%D=a(ones(nT,1),ones(maxs,1)); % Distance to climate stations
%W=a(ones(nT,1),ones(maxs,1)); % Which climate stations--cols of C
OB=a(ones(nT,1),ones(maxs,1)); % B-operator order for best model
OF=a(ones(nT,1),ones(maxs,1)); % F-operator order for best model
K2=a(ones(nT,1),ones(maxs,1)); % any signif params in final model?
%N = a(ones(nT,1),:);  % number of stations in search radius
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

k1=[1 2];
%PT= dec2dms2(TC,k1); % change coordinates from map long-lat to regular
%PC= dec2dms2(CC,k1);
%P=gcdist(PT,PC); % calc dist (km) from each tree-site to each station
%[D,W,N,I2]=near3(P,ds,maxs,maxs,1);
% D, W, N defined before;  I2 is pointer to valid chronologies, 
%		in other word, chronologies with at least one climate station
%		in the search radius
%
% Note: I did the above outside this function, and saved
% D W N I2 in distinfo

%*************** LOOP OVER CHRONOLOGIES ************************
pack

m = (1:mT)';

% Build a matrix of OE model structures.  These will be the
% candidate models.
nn=struc3(1:nb,0:nf,0,ntot); % no-delay models only


% Normally, analyze all mT chronologies; but for debugging and
% detailed analysis, might want to do a small subset, identified
% as a subset of cols of T
if O(1)==2,
	m=input('row vector telling which chronologies to do')
	if any(m<1) | any(m>mT),
		error('m out of range')
	end
end

% Start the loop over chronologies
for n =1:nT; %1:nT;
	ns=N(n); % number of stations found for this chronology
	C1=C(:,W(n,1:ns)); % col sub-matrix of the "near" stations
	z= T(:,n); % all rows of data for this chronology
	
	% Loop over climate stations; or do just one station
	if O(2)==1,
		ndo = ns;
	else
		ndo=1;
	end
	
	for k = 1:ndo;
		disp(['On chron ',int2str(n),' station ',int2str(k)])
		if O(2)==1;
			c = C1(:,k);
		else
			c= C(:,w(n));
		end
		% Get the valid data for common period; need input climate
		% series u and output chronology y, and year vector yrs
		[u,y,yrs]=overlap(c,z,flipud(YRS));
		% Get correlation coefficient between u and y
		rr1= corrcoef([y u]);
		r1=rr1(1,2);
		nyrs = length(u); % number of years for OE modeling
		%  OE-model the series with split-sample validation
		Z=[y u];
		[s,r2,ob,of,f1,f2,way,k2]=oesplit(Z,yrs,nn,nranked,O(1));
		% O(1)==1 means speedy mode; ==2 means diagnostic mode

		% Store the results for this chronology/station
		K2(n,k)=k2;
		OB(n,k)=ob;
		OF(n,k)=of;
		G(n,k)=f1;
		H(n,k)=f2;
		R1(n,k)=r1;
		R2(n,k)=r2;
		S(n,k)=s;
		NY(n,k)=nyrs;
		Q(n,k)=way;
	end

end; % of "n=" loop over chronologies


set2=' K2 OB OF R1 R2 G H S NY ';

% Store result in a .mat file of the user's choice
disp('Example file name for next prompt is sop1104a.mat')
disp(' meaning Screen.m Output Pcp, Nov-April, station subset A')
disp(' press return')
pause
[file1, newpath] = uiputfile('*.mat', 'Save Results As');
eval(['save ',newpath,file1,set2])
