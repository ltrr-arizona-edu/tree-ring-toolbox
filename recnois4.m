function Y=recnois4(yc,E,yrsc,nsim,kopt)
% recnois4:  add bootstrap or other noise to calib-period part of reconstruction
% CALL: Y=recnois4(yc,E,yrsc,nsim,kopt);
%
%************  IN ***********
%
% yc (m1 x 2)r calibration period reconstructed data, year in col 1
% E (nyrc x 1)r  calibration-period noise computed from cross-validation estimates
% yrsc (1 x 2)i  first, last years of calibration period, which also are
%    the first and last years for E
% nsim (1 x 1)i number of simulation series to generate
% kopt(1 x 2)i options
%		kopt(1) how to build the noise
%			==1 bootstrap 
%			==2 random normal, using mean and variance of calib noise
%     kopt(2) backtransform noise-added series  from log to original units 
%        ==1  no
%        ==2  yes
%
%*******   OUT ********************
% User prompted for file to store the following matrix:
%   Y (mY x nY)r  the reconstruction plus noise.  Each column a simulation
%		
%********** NOTE ***************
%
% User typically would have run recflow1.m before this.  
%
% User is given screen prompt for option to un-logtransform the final noise 
% series.  This useful if reconstruction was in terms of log10 data
%
% Noise-added series covers only the calibration period


str='NOISE INFORMATION';

%--------------  DO YOU WANT RECON+NOISE SERIES TO BE BACKTRANSFORMED 
% FROM LOG10 TO ORIGINAL UNITS
%backt = questdlg('Backtransform Noise series from log10 to original units?');
% commented out the above 11-23-99.  modified to specify as kopt(2)


%-------------- CHECK INPUT

%yo and ym
[mtemp,ntemp]=size(yc);
if ntemp~=2;
   error('yc must be 2-col tsm with year in col 1');
end

% calibration period
[mtemp,ntemp]=size(yrsc);
if mtemp~=1 | ntemp~=2,
   error('yrsc must be 1 x 2');
end
yrc = (yrsc(1): yrsc(2))'; % year vector for calibration period
nyrc = length(yrc);


mtemp=length(E);
if mtemp~=nyrc;
   error('E row size must equal nyrc');
end
% nyrc is number of calibration period years
% yrc is calibration period year vector

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
if mtemp~=1  | ntemp~=2;
   error('kopt must be 1 x 2');
end
if kopt(1) <1 | kopt(1)>2;
   error('kopt(1) must be 1 or 2');
end
if kopt(1)==1;
   str=char(str,'Noise computed by bootstrap method');
elseif kopt(1)==2;
   str=char(str,'Noise computed by sampling random normal distribution');
end
if kopt(2) <1 | kopt(2)>2;
   error('kopt(2) must be 1 or 2');
end
if kopt(2)==1;
   str=char(str,'No backtransforming of noise-added series');
elseif kopt(2)==2;
   str=char(str,'Noise-added series backtransformed from log10 units');
end
if kopt(2)==2;
   backt='Yes';
else;
   backt='No';
end;


%----------  ALLOCATE TO STORE NOISE ADDED SERIES
 % initialize Y as the calibration period reconstructon
Y = repmat(yc(:,2),1,nsim);


e=E; % cross-validation compute model noise
nyr = nyrc; % will need to generate this many years of noise
ngen=nyr;

%-------- Generate noise
if kopt(1)==1; % want bootstrap method
   % First need to stack calib pd noise to longer sequence than recon segment
   [F,I]=bootstrp(nsim,'nullfun',e);
   F=F';
   F = F(1:ngen,:);
   clear I
elseif kopt(1)==2; % random normal
   stde=std(e);
   F = normrnd(0,stde,ngen,nsim);
end
%************  Add noise to recon
Y= Y + F;

%***********************  POSSIBLY BACK TRANSFORM

if strcmp(backt,'Yes');
   Y = exp(log(10)*Y);
else
end
