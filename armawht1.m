function [Y,mods,varaut,Q,A3,SE2]=armawht1(X,yr,YRS,maxparam,PQ,kopt)
% armawht1: residuals from ARMA modeling of time series in matrix
% CALL: [Y,mods,varaut,Q,A3,SE2]=armawht1(X,yr,YRS,maxparam,PQ,kopt);
%
% Meko 3-18-98
%************* IN
%
% X (nyr x nser)r  time series in a matrix, nyr years and nser series
% yr (nyr x1)i year vector for X
% YRS (nser x 2)i  start and end year of desired period for estimating the model
%   for each series.  If [], model is fit to the full non-NaN part of each series
% maxparam (1 x 1)i  maximum number of parameters to allow in AR or ARMA model
%   Used only if kopt(3)==1
% PQ :  orders p and q of AR and MA component for each series
%   If kopt(3)==1, PQ is [], and FPE selects model structure and order
%   If kopt(3)==2, PQ is a row vector specifying p and q for all series (same model)
%   if kopt(3)==3, PQ is (nser x 2)i, with possibly different secified model for each series
% kopt (1 x 3)  options
%   kopt(1) period you want residuals for
%      ==1 Y to have residuals for all years possible given data coverage in X
%          (allows for fitting model to modern period and filtering earlier data
%           by that model
%      ==2 Y to have residuals for only the fit period specified in YRS.  All other
%         elements of Y set to NaN
%   kopt(2): model structure to consider in fitting data by FPE(used only if kopt(3)==1)
%      ==1 AR and ARMA
%      ==2 AR only
%   kopt(3): how to pick model order
%      ==1 lowest Akaike FPE 
%      ==2 same model order for all series (overrides kopt(2) and maxparam)
%      ==3 model order specified for all series (overides kopt(2) and maxparam)
% 
%
%************ OUT
%
% Y (nyr x nser)r  the whitened version of series in X
% mods{1x nser}s   the model structures (e.g., {'ARMA(1,1)','AR(2)' ,...})
% Q (nser x 2)r  portmanteau stat and its p-value for each series
% A3 (nser x 3)r  lags 1-3 of acf of model residuals
% SE2 (nser x3)r  two times the large-lag standard error for A3 
%
%
%************ NOTES
%
% Null model.  A series may have so little modelable persistence that
% ARMA modeling is not justified.  This is checked by Akaike's FPE.
% If the Akaike's FPE for all models is greater than the original variance
% of the series (computed with N-1 in the denominator), the series is 
% accepted as random already, and the original series is returned as the
% whitened series.  If so, other output assumes these values:
%  -mods   "null"
%  -varaut 0
%  -Q, A3, SE2   computed on original series instead of AR or ARMA residuals
% 
% Leading years.  One or more leading years of the time series will be lost in
% whitening, depending on order of the model.  This function changes those startup
% values to NaN based on order of the AR term.  If AR(1), one value is lost, etc.
%
%******************** CHECK INPUT ********************

%----------- X, YRS, yr
[nyr,nser]=size(X);
[m1,n1]=size(YRS);
[m2,n2]=size(yr);
if n2~=1 | m2~=nyr;
   error('yr must be cv same row size as X');
end

if isempty(YRS);
   % no problem; just means will fit models to all years of data in X
else
   if nser~=m1;
      error('row size of YRS must equal col size of X');
   end
   if min(diff(YRS'))<20;
      error('Need at least 20 years for arma modeling');
   end
   if  any(YRS(:,1)<min(yr));
      error('A start year in YRS is earlier than first year in yr');
   end
   if  any(YRS(:,1)>max(yr));
      error('An endt year in YRS is later than last year in yr');
   end

end

%------- maxparams
[m1,n1]=size(maxparam);
if ~(m1==1 &  n1==1);
   error('maxparam must be scalar');
end


%--------------  PQ, kopt
[m1,n1]=size(PQ);
[m2,n2]=size(kopt);
if ~(m2==1 & n2==3);
   error ('kopt must be 1 x 3');
end
if ~any(kopt(1)==[1 2]);
   error('kopt(1) must be 1 or 2');
end
if ~any(kopt(2)==[1 2]);
   error('kopt(2) must be 1 or 2');
end
if ~any(kopt(3)==[1 2 3]);
   error('kopt(3) must be 1,2 or 3');
end
if kopt(3)==1;  % FPE to select model
   if ~(m1==0 & n1==0);
      error('PQ must be [] if kopt(3)==1');
   end
elseif kopt(3)==2; % same model for all series
   if ~(m1==1 & n1==2);
      error('PQ must be 1 x 2 if kopt(3)==2');
   end
elseif kopt(3)==3; % specify model for each series
   if ~(m1==nser & n2==2);
      error('PQ must be nser x 2 if kopt(3)==3');
   end
end


%***************  ALLOCATE FOR OUTPUT

Y = repmat(NaN,nyr,nser);
Q = repmat(NaN,nser,2);
A3=repmat(NaN,nser,3);
SE2 = repmat(NaN,nser,3);
      

%*************  IF USING FPE, COMPUTE ORDERS OF CANDIDATE MODELS 
if kopt(2)==2; % AR modeling only
   NNar = [1:maxparam];
   NNarma=[];
elseif kopt(2)==1; % allow AR and ARMA
   NNar=[1:maxparam];
   k=0;
   for ma=1:(maxparam-1);
      for mc = 1:(maxparam-ma);
         k=k+1;
         NNarma(k,:) =[ma mc];
      end
   end
end
% nn and nn2 now hold the orders of the candidate ar and arma models

      
%***************** FIT THE MODELS

yr1=yr;
for n = 1:nser; % loop over time series
   
   x1 = X(:,n); % full period of x1, with year vector yr1
   
   % Get the subset of all non-NaN data 
   L1 = ~isnan(x1);
   x2=(x1(L1));
   yr2 = yr1(L1);
   
   % Get subset of non-NaN data to be used in modeling
   if isempty(YRS); % use all available NaN data
      x3=x2;
      yr3=yr2;
      if length(yr3)<20;
         error(['Series ' n ' has fewer than 20 yr data for modeling ']);
      end
   else; % use only the (non-NaN) data in YRS(n,:)
      yrs = YRS(n,:);
      L1 = yr2>=yrs(1) & yr2<=yrs(2);
      yr3=yr2(L1);
      x3=x2(L1);
      if length(yr3)<20;
         error(['Series ' n 'has fewer than 20 yr non-NaN data in specified period']);
      end
   end
   % Pointer to rows of X with any NaN data
   LX2 = yr1>=min(yr2) & yr1<=max(yr2);
   nyr2=length(yr2);
   % Get pointer to storage of model residuals in Y
   LY3 = yr1>=min(yr3) & yr1<=max(yr3);
   nyr3=length(yr3);
   
   
   % Now have:
   % x1,yr1  -- the full col of data for the series, maybe with NaNs
   % x2,yr2  -- the non-NaN part
   % x3,yr3  -- the segment to be used to fit the model
   % Note that the mean has not yet been subtracted
   
   % number of candidate models of each structure
   nar = length(NNar);
   narma = size(NNarma,1);

   % Pre allocate for FPE
   Far = repmat(NaN,nar,1);
   if narma>0;
      Farma = repmat(NaN,narma,1);
   else
      Farma=[];
   end
   
   
   u=x3; % put the series in col vector u
   
   % compute mean and standard deviation for modeling period
   mnu=mean(u);
   stdu = std(u);
      
   % ready the chron for arma modeling by subtracting mean
   u = u -mnu;
   
   
   %-----------  FIT MODELS USING FPE
   
   if kopt(3)==1; % use FPE to get best model
      % Compute and store FPE for candidate models
      for n1 = 1:nar;
         thar = ar(u,NNar(n1));
         Far(n1) = thar(2,1); % store FPE
      end
      if narma>0;
         for n2 = 1:narma;
            tharma = armax(u,NNarma(n2,:));
            Farma(n2) = tharma(2,1); % store FPE
         end
      else
      end
      
      
      % Find best ar model and its FPE
      [Fbest,i] = sort(Far);
      nnar = NNar(i(1));
      Far1 = Fbest(1);
         
      % Find best arma model
      if narma>0;
         [Fbest,i]= sort(Farma);
         nnarma = NNarma(i(1),:);
         Farma1 = Fbest(1); 
      end
   
   
      %-------- Find best overall model
      if narma>0;% If ar and arma allowed
         [F1,ii]=sort([Far1 Farma1]);
         if ii(1)==1; % ar
            modtype='AR';
            nn  = nnar; 
            strmoda = sprintf('(%1.0d)',nn);
            strmod = [modtype strmoda]; % string for model, like AR(1)
         else
            modtype='ARMA';
            nn = nnarma;
            strmoda = sprintf('(%1.0d,%1.0d)',nn);
            strmod = [modtype strmoda]; % string for model, like ARMA(1,1)
         end
      else; % only ar models allowed   
         modtype='AR';
         nn  = nnar; 
         strmoda = sprintf('(%1.0d)',nn);
         strmod = [modtype strmoda]; % string for model, like AR(1)
      end
      
   else  % fit same model for all series, or specified model
      if kopt(3)==2; %  same model all series
         p=PQ(1); q=PQ(2);  % ar and ma orders
      else
         p=PQ(n,1); q=PQ(n,2);
      end
      
      if q==0;
         modtype='AR';
         nn=p;
         strmoda = sprintf('(%1.0d)',nn);
         strmod = [modtype strmoda]; % string for model, like AR(1)
      else
         modtype='ARMA';
         nn=[p q];
         strmoda = sprintf('(%1.0d)',nn);
         strmod = [modtype strmoda]; % string for model, like AR(1)
      end
   end
   
         
   %----------- Refit model and get model residuals
   
   if strcmp(modtype,'AR');
      th=ar(u,nn);
   elseif strcmp(modtype,'ARMA');
      th=armax(u,nn);
   end
   figure(10);
   e=resid(u,th);
   close(10);
      
   % pct variance due autocorrelation
   vare = th(1,1); % estimated noise variance
   varu = var(u);  % variance of chron
   pct= 1 - vare/varu;
   pct=pct*100;
               
   % Portmanteau statistic
   [ree,Stderr2,r95]=acf(e,20); % acf of model residuals, in ree
   if strcmp(modtype,'AR');
      p = nn;
      q=0;
   else
      p = nn(1); q=nn(2);
   end
   [P,pval]=portmant(ree,nyr3,p,q,20);
   %strport{n} = sprintf('Portmanteau stat = %8.4f, P-value=%7.5f',P,pval);
   
   
   %---------- Rename some output
   mods{n}=strmod; % string identifier for model order
   Q(n,:) = [P pval]; % portmantea and p-value
   A3(n,:)=ree(1:3); % acf of residuals at lags 1-3
   varaut(n)=pct; % percent variance due to modeled autocorrelation
   SE2(n,:) = Stderr2(1:3);
   
     
   % Put residuals (with original mean restored) in Y.  
   if kopt(1)==2;  % residuals from fit period only are to be used
      Y(LY3,n)=e+mnu;
   elseif kopt(1)==1;  % arma-filtered values from outside fit period also wanted
      u2=x2-mnu; % original data (all non-NaN data) with mean for model period subtracted
      u2p=predict(u2,th,1); % prediction based on fitted arma model
      e2 = u2-u2p; % difference between observed and model predicted
      Y(LX2,n)=e2+mnu;  % add model-period mean back
   end
   
   % Replace any startup values depending on nonexistent data before the first
   % available observation with NaN
   itemp=find(LX2);  % index to rows of X with non-NaN data
   Y(itemp(1:nn(1)),n)=NaN;  % change first nn(1) values of the whitened series
      % to NaN if those values are also the first nn(1) Non-NaN data in the 
      % original series
   
   
   
   
   
   %------------ Handle null case in which FPE higher than original variance of series
   if th(2,1)>var(u);
      mods{n}='Null';
      varaut(n)=0;
      
      % Portmanteau statistic, acf
      [ru,su2,r95]=acf(u,20); % acf of original series, in ru
      [P,pval]=portmant(ru,length(u),0,0,20);
      Q(n,:)=[P pval]; % portmanteau of original series
      A3(n,:) = ru(1:3);
      SE2(n,:) = su2(1:3);
      
      % Replace residuals with original series
      if kopt(1)==2;
         Y(LY3,n) = u + mnu;
      elseif kopt(1)==1;
         Y(LX2,n) = x2;
      end
      
   end
end

   
      