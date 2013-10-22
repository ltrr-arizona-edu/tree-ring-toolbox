function Y=recnois3(yo,ym,YRS,nsim,kopt,yrsmask)
% recnois3:  add bootstrap or other noise to calibration period part of
% reconstruction to noise-added matrix for the calibration period.
% CALL: recnois3(yo,ym,YRS,nsim,kopt);
%
%************  IN ***********
%
% yo (m1 x 2)r observed time series of predictand, with year in col 1
% ym (m2 x 1)r reconstructed time series of predictand, with year in col 1
% YRS(1 x 2)i  first, last years of calib period
% nsim (1 x 1)i number of simulation series to generate
% kopt(1 x 1)i options
%		kopt(1) how to build the noise
%			==1 bootstrap
%			==2 random normal, using mean and variance of calib noise
% yrsmask(1 x ?)i  years to mask out in computing noise
%   yrsmask==[] if no years masked
%
%*******   OUT ********************
% User prompted for file to store the following matrix:
%   Y (mY x nY)r  the reconstruction plus noise.  Each column a 
%   noise-added series
%		
%********** NOTE ***************
%
% User is given screen prompt for option to un-logtransform the final noise 
% series.  This useful if reconstruction was in terms of log10 data


%-------------- CHECK INPUT

%yo and ym
[mtemp,ntemp]=size(yo);
if ntemp~=2;
   error('yo must be 2-col tsm with year in col 1');
end
[mtemp,ntemp]=size(ym);
if ntemp~=2;
   error('ym must be 2-col tsm with year in col 1');
end

% nsim
[mtemp,ntemp]=size(nsim);
if mtemp~=1 | ntemp~=1;
   error('nsim must be scalar');
end
if nsim<1;
   error('nsim must be Positive');
end

% kopt
[mtemp,ntemp]=size(kopt);
if mtemp~=1  | ntemp~=1;
   error('kopt must be 1 x 1');
end
if kopt(1) <1 | kopt(1)>2;
   error('kopt(1) must be 1 or 2');
end


%----- YRS
[mtemp,ntemp]=size(YRS);
if mtemp~=1 & ntemp~=2; 
	error('YRS must be 1 x 2');
end

% make year vectors
yr1 = yo(:,1);
yr2 = ym(:,1); 
if ~all(diff(yr1)==1)
   error('year col of yo does not increment by 1');
end
if ~all(diff(yr2)==1)
   error('year col of ym does not increment by 1');
end
yrs1=[min(yr1) max(yr1)];
yrs2=[min(yr2) max(yr2)];

yr3 = YRS(1,1):YRS(1,2); % year vector for the calibration period

% Noise period must have both observed and model recon data
if YRS(1,1) < yrs1(1) ;
   error('YRS(1,1) settings inconsistent with start years of yo or ym');
end
if YRS(1,2) > yrs1(2) ;
   error('YRS(1,2) settings inconsistent with end years of yo or ym');
end


nyr1=length(yr1);
nyr2=length(yr2);
nyr3=length(yr3);

%--------------  DO YOU WANT RECON+NOISE SERIES TO BE BACKTRANSFORMED 
% FROM LOG10 TO ORIGINAL UNITS
backt = questdlg('Backtransform Noise series from log10 to original units?');

%----------  ALLOCATE
Y = repmat(NaN,nyr3,nsim);


%----------- BUILD POINTERS

% Compute Lomask == logical pointer to years in yo to not mask in computing noise
if ~isempty(yrsmask);
   nmask = length(yrsmask);
   YRSmask = repmat(yrsmask,nyr1,1);
   YR1 = repmat(yr1,1,nmask);
   Lmask = YR1 == YRSmask;
   Ltemp = (any(Lmask'))';
   Lomask = ~Ltemp;  
else
   Lomask = ones(nyr1,1);
   Lomask=logical(Lomask);
end

% Compute Lmmask == logical pointer to years in ym to not mask in computing noise
if ~isempty(yrsmask);
   nmask = length(yrsmask);
   YRSmask = repmat(yrsmask,nyr2,1);
   YR2 = repmat(yr2,1,nmask);
   Lmask = YR2 == YRSmask;
   Ltemp = (any(Lmask'))';
   Lmmask = ~Ltemp;  
else
   Lmmask = ones(nyr2,1);
   Lmmask=logical(Lmmask);
end



Lco= yr1>=YRS(1,1) & yr1<=YRS(1,2);  % to calibration period in yo
Lcon = Lco & Lomask; % to calibration years in yo actually to be used for noise calc
Lcm=yr2>=YRS(1,1) & yr2<=YRS(1,2);  % to calibration period in ym
Lcmn = Lcm & Lmmask; % to calib perd yrs in ym actually to be used for noise calc
ngen = length(yr3);


%*******************  Compute the calibration period noise series
yoc = yo(Lco,2); yocn = yo(Lcon,2); % observed, with all data, and with outliers
   % ommitted, depending on yrsmask
ymc = ym(Lcm,2); ymcn = ym(Lcmn,2); % likewise for reconstructed data in calib pd
ec = yocn-ymcn;  % calibration period noise: observed minus model
nc = length(ec);

%************* Screen check on explained variance
clc;
ev = 1 - var(ec)/var(yocn);
disp('Quality control: check the following variance explained fraction'); 
disp('against your regression R-squared');
if isempty(yrsmask);
   disp('Noise was computed using all calibration years');
else
   disp('Noise was computed after omitting the following years');
   disp(yrsmask);
end
strscr1 = sprintf('EV = %5.2f\n',ev);
disp(strscr1);
if kopt(1)==1;
   disp('Noise computed by bootstrap method');
elseif kopt(1)==2;
   disp('Noise computed by sampling random normal distribution');
end


%*****************  Pull noise-free recon segment 
x = ym(Lcm,2); % this is the recon data that will add noise to


%*********** COMPUTE NOISE SERIES FOR CALIBRATION PERIOD

if kopt(1)==1; % want bootstrap method
   % First need to stack calib pd noise to longer sequence than recon
	nstack = ceil (ngen/nc);  
	e3 = repmat(ec,nstack,1);
	[E,I]=bootstrp(nsim,'nullfun',e3);
	E=E';
	E = E(1:ngen,:);
	clear I
elseif kopt(1)==2; % random normal
	stde=std(ec);
	E = normrnd(0,std(ec),ngen,nsim);
end


%************  ADD CALIBRATION PERIOD NOISE TO CALIB-PERIOD RECONSTRUCTION
X = repmat(x,1,nsim);
X = X + E;
Y = X; % copy-rename  recon+ noise matrix

%***********************  POSSIBLY BACK TRANSFORM
if strcmp(backt,'Yes');
   Y = exp(log(10)*Y);
else
end


%************************ STORE NOISE
[file1,path1]=uiputfile('noical*.mat','Store noise here');
pf1=[path1 file1];
eval(['save ' pf1 ' Y;']);


