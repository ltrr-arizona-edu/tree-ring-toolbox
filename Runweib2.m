function runweib2
% runweib2: Index drought summary of reconstruction simulations vs actual 
% CALL: runweib2;
%
% Written for Jordan precipitation reconstruction. 
% Method follows Sadeghipour, J. and Dracup, J. A., 1985,  'Regional frequency 
% analysis of hydrologic multiyear droughts', Water Resources Bulletin,
% 21(3), 481-487. runweib2.m is special request from ramzi because needed figure
% that would be clear when greatly reduced ( 11/18/98).  Otherwise, runweib2 is same
% as runweib1
% 


%***************** STEPS **************************************
%
% Get .mat file with runsnoi1.m output summary runs data
% Loop over data types: actual instr pd, recon instr pd, long-term recon
%		Convert the drought severities to standardized drought severities
%		Sort the events from most severe to least
%		Compute the Weibull probabitie; store these and the sorted stdzd
% End Loop
%
% Convert the Simulation severities to standardized severities
% Count number of events (all lengths) in each simulation; compute median 
% Compute Weibull probabilities, using the median number of events
% Compute and store boxplot info for each Weibull prob
%
% PLOTS...
%
% x-axis: excedance prob
% y-axis: standardized severity
% Points: (1) dots for actual data (2) boxplots for noisy reconstruction

%--------------  GET FILE WITH RUNSNOI1.M OUTPUT

[file1,path1]=uigetfile('*out.mat','Input data from a runsnoi1.m output');
pf1=[path1 file1];
eval(['load ' pf1]);


%-----------  COMPUTE STANDARDIZED SERIES FOR ALL EXCEPT SIMULATIONS

% Note that computations for all series use mean and std dev from the
% actual data

z1 = (s1 - d1 * mn1) ./ (sqrt(d1 * var1)); 
z2 = (s2 - d2 * mn2) ./ (sqrt(d2 * var2)); 
z4 = (s4 - d4 * mn4) ./ (sqrt(d4 * var4)); 

% Sorted (descending severity) and Weibull probs

z1sort = flipud(sort(z1));
num1 = length(z1);
rankvect = (1:num1)';
weib1 = rankvect / (num1+1);

z2sort = flipud(sort(z2));
num2 = length(z2);
rankvect = (1:num2)';
weib2 = rankvect / (num2+1);

z4sort = flipud(sort(z4));
num4 = length(z4);
rankvect = (1:num4)';
weib4 = rankvect / (num4+1);


%------------  STANDARDIZE THE SIMULATION SEVERITIES
[mS,nS]=size(S);

Z = (S -  D .* (repmat(mnns,mS,1)))  ./ sqrt(D .* repmat(varns,mS,1));

% Sort
Zs=sort(-Z);
Zs=-Zs;

% Get median number of events in the noisy series
Ltemp = ~isnan(Zs);
numev = median(sum(Ltemp));
nummax = max(sum(Ltemp));

W= Zs(1:numev,:) ; % truncate 
W=W';

% Get median standardized drought for each rank
wmed = median(W);

% Compute Weibull 
p = (1:numev)' / (numev);


% Compute Z values for the historical events held in d1 and s1
Zinst = (s1 - d1*median(mnns)) ./ sqrt(d1 * median(varns));

pinst = interp1(wmed',p,Zinst); % probability points



% Plot the Exceedance probabilities of standardized severity
figure(1) 
h1 = plot(p,wmed,pinst,Zinst,'o');
grid
xlabel('Exceedance Probability');
ylabel('Standardized Severity');
title('Exceedance Probability of Drought Severity');




%***************  PLOTS OF SEVERITIES OF MOST SEVERE DROUGHT

figure(2)
hold off

maxj=max([d1;d4]); % duration of longest reconstructed drought in noise-free recon
%       and observed ppt record

% Initialize storage for median, .05 quantile and .95 quantile of severity for the most
% severe drought of duration 1-yr, 2-yr, etc.  Note that there will be 1000
% "most severe" n-year droughts, one for each simulation, and that the quantiles
% will be computed from those 1000 values
Smed=repmat(NaN,1,maxj);
LO=repmat(NaN,1,maxj);
HI=repmat(NaN,1,maxj);

% Initialize storage for number of simulations (of a possible 1000) that have
% at least one n-yr drought, for n = 1 yr, 2 yr, etc.
Ntot = zeros(maxj,1);

highx=0;
longest = nanmax(nanmax(D));  % longest run in any simulation

k1=1; % while control
j=0; % counter for the drought current drought duration 
while k1==1;
   j=j+1; % duration (yr)
   
   % First we handle the noise-added recons
   Stemp = S;  % copy-store matrix of severities
   Ltemp1 = D==j; % pointer to elements of D (and of Stemp) for the j-yr droughts
   Stemp(~Ltemp1)=NaN; % replace elements of severity matrix with NaN for all
   %   except the duration of interest
   ntot = sum(Ltemp1); % row vector of number of j-yr droughts in each simulation
   Ntot(j)=sum(ntot>0); %  number of simulations with at least one j-yr drought
   fnone = sum(ntot==0)/nS; % decimal fraction of sims with no drought of length j 
      
   if fnone >0.95 | j>longest;
      % If more than 95% of the simulations fail to yield any runs with
      % this run length, or if j exceeds the longest duration in any
      % simulation
      k1=0;
      break
   else; % at least 5% of the simulations have a drought of duration j
      highx = j; % copy-store duration
      Smax = nanmax(Stemp); % row vector of maximum severity of any j-yr drought
      %   in each simulation
      Lkeep = ~isnan(Smax); % pointer to simulations with at least one j-yr drought
      Smax = Smax(Lkeep); % row vector of maximum severities, with NaNs gone
      xbars = prctile(Smax,[5 50  95]); % percentiles of maximum severity
      Smed(j) = xbars(2);  % the median value of the maximum drought severity
      % Compute quantities needed for use of errorbar.m
      LO(j)=Smed(j)-xbars(1); % difference between that median and the 5th percentile
      HI(j)=xbars(3)-Smed(j); % diff between the 95th percentile and the median
   end
   
   % Next we handle the instrumental record
   Stemp=s1; % copy-store the col vector of severities
   Ltemp1 = d1==j; % point to rows with droughts of duration j
   if ~any(Ltemp1); % if no droughts of this duration
      skey(j)=NaN;
   else; % if at least one j-yr drought
      Stemp(~Ltemp1)=NaN; % substitute NaN for severity for droughts of other duration
      Smax = nanmax(Stemp); % maximum severity of j-yr drought
      skey(j) = Smax; % store maximum severity for this duration, j yr
   end
   
   
end

%------ Start the figure
% errorbars from 5th to 95th percentile, with line graph thru medians
hpe=errorbar((1:highx),Smed(1:highx),LO(1:highx),HI(1:highx));
hold on
hp2a=plot((1:highx)',skey(1:highx),'.','MarkerSize',32); % max severity for observed data
hp2b=plot((1:highx),Smed(1:highx),'o','MarkerSize',12); % median of max severity for
%       the 1000 simulations
title1='Most Severe Droughts for Various Run Lengths';
title2=[int2str(nS) ' Simulations']
title({title1,title2});
xlabel('Duration (yr)','FontSize',16);
ylabel('Severity (mm )','FontSize',16);
set(gca,'Xtick',[1:highx]);

set(hpe(1),'LineWidth',3);
set(hpe(2),'LineWidth',2);
set(hp2b,'LineWidth',2);

% Anotate with number of simulations having any runs of this length
for k = 1:highx;
   txt1 = int2str(Ntot(k));
   xpoint = k;
   ypoint = HI(k)+Smed(k);
   if Ntot(k)<nS;
      htext=text(xpoint,ypoint,txt1);
      set(htext,'VerticalAlignment','bottom',...
         'HorizontalAlignment','center',...
         'FontSize',14)
   end
 end
 
 Xlim2=get(gca,'Xlim');
 Ylim2=get(gca,'Ylim');
 
 set(gca,'YTick',[50:100:550],'Fontsize',14);
 
 
grid


%************ NOISE-FREE RUNS SUMMARY **********

hold off
figure(3);


D4 = repmat(d4,1,highx);
itemp = 1:highx;
I4 = repmat(itemp,length(d4),1);
S4 = repmat(s4,1,highx);
L4 = D4==I4;
S4(~L4)=NaN;
maxrun4 = nanmax(S4);

D2 = repmat(d2,1,highx);
itemp = 1:highx;
I2 = repmat(itemp,length(d2),1);
S2 = repmat(s2,1,highx);
L2 = D2==I2;
S2(~L2)=NaN;
maxrun2 = nanmax(S2);


plot((1:highx)',skey(1:highx),'.','MarkerSize',22);
xlabel('Run Length (yr)');
ylabel('Run Sum (mm )');
set(gca,'Xtick',[1:highx]);
hold on
plot((1:highx),maxrun4,'o','Markersize',12);
line((1:highx),maxrun4);

plot((1:highx),maxrun2,'square','Markersize',14);

set(gca,'Xlim',Xlim2,...
   'Ylim',Ylim2);
title('Highest Run Sums of Precipitation Deficit');
legtext = [...
   'Observed, 1946-95        ';...
   'Reconstruction, 1946-95  ';...
   'Reconstruction, 1600-1995'];
h = get(gca,'Children');
hleg = legend([h(4);h(1);h(3)],legtext);
grid

