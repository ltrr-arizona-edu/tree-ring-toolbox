function Y=recnois2(yo,ym,E,yrsc,YRSr,nsim,kopt)
% recnois2:  add bootstrap or other noise to sub-period reconstructions
% CALL: Y=recnois2(yo,ym,E,yrsc,YRSr,nsim,kopt);
%
%************  IN ***********
%
% yo (m1 x 2)r observed time series of predictand, year in col 1
% ym (m2 x 2)r reconstructed time series of predictand, year in col 1
% E (nyrc x nper)r  or (1 x 3)r 
%   if kopt(1)==1 calibration-period noise computed from cross-validation estimates (nyrc x nper)
%   if kopt(1)==2 RMSE of crossvalidation (1 x nper)
% yrsc (1 x 2)i  first, last years of calibration period, which also are
%    the first and last years for E
% YRSr(nper x 2)i  first, last years of reconstruction subperiods
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
% Assumes same calibration period used for all sub-period models
% Assumes ym covers only the reconstruction periods, and covers all
%   the reconstruction sub-periods listed in YRSr
%
% Observed predictand data (without noise) is tacked onto the tail end
% of the noise+recon data.  Thus user will see that last few rows of Y will
% all have identical data
%
% User is given screen prompt for option to un-logtransform the final noise 
% series.  This useful if reconstruction was in terms of log10 data
%
% Noise-added series covers entire length of reconstruction and observed data


str='NOISE INFORMATION';

%--------------  DO YOU WANT RECON+NOISE SERIES TO BE BACKTRANSFORMED 
% FROM LOG10 TO ORIGINAL UNITS
%backt = questdlg('Backtransform Noise series from log10 to original units?');
% commented out the above 11-23-99.  modified to specify as kopt(2)


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
nyrm=size(ym,1);

% calibration period
[mtemp,ntemp]=size(yrsc);
if mtemp~=1 | ntemp~=2,
   error('yrsc must be 1 x 2');
end
yrc = (yrsc(1): yrsc(2))'; % year vector for calibration period
nyrc = length(yrc);


% E
[mtemp,nper]=size(E);
if kopt(1)==1; % bootstrap: need actual error values
    if mtemp~=nyrc;
        error('E row size must equal nyrc if using bootstrap');
    end
elseif kopt(1)==2; % monte carlo
    if mtemp~=1;
        error('E row size must be 1 if using Monte Carlo');
    end;
end
% nper is the number of sub-periods
% nyrc is number of calibration period years
% yrc is calibration period year vector


% YRSr 
[mtemp,ntemp]=size(YRSr);
if ntemp~=2 | mtemp~=nper;
   error('YRSr row size must equal nper, and col size must be 2');
end
if YRSr(1,1)~=min(ym(:,1)) | YRSr(nper,2)~=max(ym(:,1));
   error('Mismatch in first or last year of ym and YRSr');
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
yrs4=[min(ym(:,1))  max(yo(:,1))];  % from earliest recon year to last observed
yr4=(yrs4(1):yrs4(2))';
nyr4 = length(yr4);
Y = repmat(NaN,nyr4,nsim);


%--------- Fill first part of Y with the noise-free reconstruction, last part with observed
%  last part with observed
Y(1:nyrm,:)=  repmat(ym(:,2),1,nsim);
Y((nyrm+1):nyr4,:)= repmat(yo(:,2),1,nsim);



%********************* LOOP OVER SUB-PERIODS

for n = 1:nper;
   yrs=YRSr(n,:);
   yr = (yrs(1):yrs(2))';
   nyr = length(yr); % will need to generate this many years of noise
   ngen=nyr;
   L1 = yr4>=yrs(1) & yr4<=yrs(2); % pointer to rows in Y that noise will go
   e = E(:,n);  % calibration period noise series for this sub-period model, if kopt(1)==1; or 
   %  RMSE of cross-validation (a scalar) for the subperiod model (if kopt==2);
   
   %-------- Generate noise
   if kopt(1)==1; % want bootstrap method
      % First need to stack calib pd noise to longer sequence than recon segment
      nstack = ceil (ngen/nyrc);  
      e3 = repmat(e,nstack,1);
      [F,I]=bootstrp(nsim,'nullfun',e3);
      F=F';
      F = F(1:ngen,:);
      clear I
   elseif kopt(1)==2; % random normal
      stde=e; % new version
      %stde=std(e); J. AWRA VERSION
      F = normrnd(0,stde,ngen,nsim);
   end
   
   %************  Add noise to recon
   Y(L1,:)= Y(L1,:) + F;
end

%***********************  POSSIBLY BACK TRANSFORM

if strcmp(backt,'Yes');
   Y = exp(log(10)*Y);
else
end


%************************ STORE NOISE
% 11-23-99 I commented next 3 lines out & did storing within a calling script
%[file1,path1]=uiputfile('nois*.mat','Store noise here');
%pf1=[path1 file1];
%eval(['save ' pf1 ' Y;']);


