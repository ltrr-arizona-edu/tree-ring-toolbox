function runs01(pfs,jcol,kopt,YRS,colors,textinf,zcrit)
% runs01:  runs analyis and runs plot summary for a time series   
% CALL: runs01(pfs,jcol,kopt,YRS,colors,textinf,zcrit);
%
% Meko 4-28-97
%
%****************** IN ***********************************
%
% pfs (1 x ?)?  path & filename of .mat file holding time series to be summarized
%   Example: 'c:\projs\ai1\recpcp.mat'
% jcol (1 x 1)i;  which column in the files pfs the data should come from
%   Example: 2 means get data col 2 of file named by pfs
% kopt(1 x 1)i  options
%   kopt(1) ==1 departure is defined as value below threshold
%           ==2 departure is defined as value above threshold
%   kopt(2) ==1 plot anomalies in original units
%           ==2 plot standardized anomalies
% YRS (1 x 2)i>  start, end year of:
%   row 1: years of data to produce bar chart for -- can be a subset of years of full data
%   row 2: years to compute monthly means and std devs on
% colors(2 x 3)i  colors for upper and lower anomalies
% textinf {}  text information
%	{1} type of variable: 'PCP', 'TMP', or 'FLOW'
%	{2} time unit of data : e.g., 'OCT-APR'
%	{2} units:  e.g., 'mm'
%	{3} namestn: station name: e.g., 'JORDAN'
% zcrit (1 x 1)r  threshold critical value from which to measure a departure
%
%
%******************** OUT **************************************
%
% In figure 1: barcharts of the anomalies from long-term seasonal means
%   for each season, with annotated runs summary
%
%
%*********************** NOTES *************************************
%
% Horizontal line is plotted at zero, indicating departure zero

close('all'); % close figure windows if any are open

%------------ UNLOAD AND CHECK TEXT INPUT
datatype = textinf{1};
season = textinf{2};
units = textinf{3};
namestn = textinf{4};


% Check YRS
[mtemp,ntemp]=size(YRS);
if ~(mtemp==2 & ntemp==2);
   error('YRS should be 2 x 2');
end

% Check out kopt
if length(kopt)~=2;
   error('kopt should be 1 x 2');
end


a=NaN;

% Initialize  matrices to hold the data for statistical and plot periods
yrstats=(YRS(2,1):YRS(2,2))'; % year vector for stats period
nstats=length(yrstats);
yrplot=(YRS(1,1):YRS(1,2))';  % for plot period
nplot=length(yrplot);
Z1=repmat(NaN,nplot,1);
Z2=repmat(NaN,nstats,1);


%************** 

clear X Y Z

% Get the .mat file that has matrix
j=jcol-1; % the col of the matrix in pfs that will have the desired data,
% after neglecting the year col
disp([  '   Series ' char(namestn) ' from ' char(pfs) ', variable # ' int2str(j)]);
eval(['load ' char(pfs)]);
% Figure out if data file is X or Y, or Z, and make it X
if exist('X')~=1;
   if exist('Y')==1;
      X=Y;
   elseif exist('Z')==1;
      X=Z;
   else
      error('Input .mat file has no matrix by name of X,Y or Z');
   end
end
   
% Get year vector for X
[mX,nX]=size(X);
yr = X(:,1);
if nX<jcol;
   error('Too few cols in X for a year and specified jcol entry');
end
   
% Lop year col off X
X(:,1)=[];
   
% Check out YRS
yrgo1=yr(1);
yrsp1=yr(mX);
if any(YRS(:,1)<yrgo1);
   clc
   disp(char(pfs));
   error('Starting years in YRS inconsistent with years in data matrix');
end
if any(YRS(:,2)>yrsp1);
   clc
   disp(char(pfs));
   error('End years in YRS inconsistent with years in data matrix');
end
   
% Pull data for plot period 
L1= yr>=YRS(1,1) & yr<=YRS(1,2); % years to plot
z1=X(L1,j);
   
% Pull data for stats period and put in col of Z2
L2= yr>=YRS(2,1) & yr<=YRS(2,2); % years to plot
z2=X(L2,j);
   
% Check that no NaN data in plot period 
if any(isnan(z1));
   error('NaN data in period for analysis and plotting');
end
   
% Get monthly 'long term' mean and std dev for stats period
xmean=nanmean(z2);
xstd = nanstd(z2);
   
% Rename cv of plot series for convenience
z3=z1;  
   
nobs=length(z3); % number of years of data in analysis period
   
%---------  If desired, convert series and critical level to standardized anomalies
if kopt(2)==2;
	z3 = (z3-xmean)/xstd;
	zcrit = (zcrit -xmean)/xstd;
	units='stdzd units';
else
end


% Form departures from critical level
z3=z3-zcrit;
   
% Compute runs statistics
if kopt(1)==1;
	ksign=-1; % interested in drought, not wetness
else
	ksign=1;
end

xc=0;  % threshold in terms of anomaly from critical level is zero

%------ COMPUTE RUNS
[Px,nx,sx]=runs(z3,yrplot,xc,ksign);
   
% Get ending year of runs
yearend = Px(:,2);
   
%******************* Start the figure
figure(1)
   
% Text info summarizing runs with largest run sum
% If a tie for largest run-sum, select that with shortest duration (most intense)
[sum1,isum1]=max(sx);
Ltemp=sx==sum1;
f1=find(Ltemp);
[temp,itemp]=min(nx(Ltemp));
isum1=f1(itemp);
sum1=sx(isum1);
long1=nx(isum1);  % run-length of run with greatest run-sum
str2a=sprintf('%g, ending in %4.0f',sum1,yearend(isum1));
str3a=sprintf(' (L = %4.0f)',long1);
   
% Text info summarizing runs with largest run length
% If a tie, select that with the greatest run-sum
[num1,inum1]=max(nx);
Ltemp=nx==num1;  % maybe have ties for longest
f1=find(Ltemp); % index to ties for longest
[temp,itemp]=max(sx(Ltemp)); 
inum1=f1(itemp);
num1=nx(inum1);
sum2=sx(inum1); % run-sum of run with longest length
str2b=sprintf('%g yr, ending in %4.0f',num1,yearend(inum1));
str3b=sprintf(' (S = %g %s)',sum2,units);
   
% Make separate series for those above zero and below zero
z3hi = z3;
z3hi(z3<=0)=NaN;
h1=bar(yrplot,z3hi);
set(h1,'FaceColor',colors(1,:));
hold on
   
z3lo =z3;
z3lo(z3>=0)=NaN;
h2=bar(yrplot,z3lo);
set(h2,'FaceColor',colors(2,:));
   
% Build y label
ylabel(['Departure (' units ')']);
xlabel('Time')
set(gca,'XGrid','on');
   
% Annotate with year
ybot=get(gca,'YLim');
ybot=ybot(1);
% for n = 1:nplot;
 %     text(n*nfiles-(nfiles-1),ybot-.20,int2str(yrplot(n)));
  % end
   
% Add info on run sum and run length
pltext(.3,.95,8,['S_{max} = ' str2a str3a]);
pltext(.3,.90,8,['L_{max} = ' str2b str3b]);
   
% Title
titlab1=sprintf(' %4.0f-%4.0f',YRS(2,1),YRS(2,2));
title([char(namestn) ' ' datatype ' Departures' ' from ' titlab1 ' Mean']);
   
hold off
