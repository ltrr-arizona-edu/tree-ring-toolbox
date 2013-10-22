function pltdep2(pfs,jcols,kopt,YRS,colors,seasons)
% pltdep2:  Seasonal bar chart of time series of departures, with runs summary 
% CALL: pltdep2(pfs,jcols,kopt,YRS,colors,seasons);
%
% Meko 4-28-97
%
%****************** IN ***********************************
%
% pfs{?)}c paths & filenames of .mat files holding seasonalized climate series
%   Example: {'d:\outfls\coolbig.mat','d:\outfls\warmbig.mat'}... for 2-season
% jcols (1 x ?)i;  which column in the files pfs the data should come from
%   Example: [2 2] means get coolbig data from col 2, warmbig data from col 2
% kopt(1 x ?)i  options
%   kopt(1) ==1 ppt
%           ==2 temperature
%   kopt(2) ==1 plot anomalies in original units
%           ==2 plot standardized anomalies
% YRS (1 x 2)i>  start, end year of:
%   row 1: years of data to produce bar chart for -- can be a subset of years of full data
%   row 2: years to compute monthly means and std devs on
% colors(2 x 3)i  colors for upper and lower anomalies
% seasons {?} contains char strings naming the seasons represented by each file in pfs
%   Example: {'NV-AP','MY-OC'}
%
%******************** OUT **************************************
%
% In figure 1, a bar chart of the anomalies from long-term seasonal means
% In figure 2, a bar chart of long-term seasonal means
%
%
%*********************** NOTES *************************************
%
% Each seasonal time series is assumed to be in a different file, as named in
%    pfs.  The col holding the seasonal series in the i-th file is given by 
%    jcols(i).
% Horizontal line is plotted at zero, indicating departure zero



% Find out how many data files will need (assumes each season's data from separate file)
nfiles=size(pfs,2);

% Check season names
ntemp=size(seasons,2);
if ntemp~=nfiles;
   error('seasons should be same dim-2 size as pfs');
end

% Check YRS
[mtemp,ntemp]=size(YRS);
if ~(mtemp==2 & ntemp==2);
   error('YRS should be 2 x 2');
end

% Check out kopt
if length(kopt)~=2;
   error('kopt should be 1 x 2');
end
if kopt(1)==1; % precip
   dtype='PPT';
   units='in';
elseif kopt(1)==2;
   dtype='TMP';
   units='^oF';
else
   error('kopt(1) must be 1 or 2');
end



% Initialize  matrices to hold the seasonalized data for the nfiles
a=NaN;
yrstats=(YRS(2,1):YRS(2,2))'; % year vector for stats period
nstats=length(yrstats);
yrplot=(YRS(1,1):YRS(1,2))';  % for plot period
nplot=length(yrplot);
Z1=a(ones(nplot,1),ones(nfiles,1));
Z2=a(ones(nstats,1),ones(nfiles,1));


%**************  LOOP OVER SEASONAL FILES
clc
disp('Looping over seasonal files to build data matrices');

for n =1:nfiles;
     clear X Y Z
   pf1=cellstr(pfs(1,n)); % path and filename
   seasnm=cellstr(seasons(1,n)); % name of season -- used in labeling
   % Get the .mat file that has matrix
   j=jcols(n)-1; % the col of the matrix in pf1 that will have the desired data,
   % after neglecting the year col
   disp([  '   Season ' char(seasnm) ' from ' char(pf1) ', variable # ' int2str(j)]);
   eval(['load ' char(pf1)]);
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
   if nX<jcols(n);
      error('Too few cols in X for a year and specified jcols entry');
   end
   
   % Lop year col off X
   X(:,1)=[];
   
   % Check out YRS
   yrgo1=yr(1);
   yrsp1=yr(mX);
   if any(YRS(:,1)<yrgo1);
      clc
      disp(char(pf1));
      error('Starting years in YRS inconsistent with years in data matrix');
   end
   if any(YRS(:,2)>yrsp1);
      clc
      disp(char(pf1));
      error('End years in YRS inconsistent with years in data matrix');
   end
   
   % Pull data for plot period and put in col of Z1
   L1= yr>=YRS(1,1) & yr<=YRS(1,2); % years to plot
   Z1(:,n)=X(L1,j);
   
   % Pull data for stats period and put in col of Z2
   L2= yr>=YRS(2,1) & yr<=YRS(2,2); % years to plot
   Z2(:,n)=X(L2,j);
   
end

% Check that no NaN data in plot period for any series
if any(any(isnan(Z1)));
   error('NaN data in period for analysis and plotting');
end


% Get monthly 'long term' means, medians, and std devs
xmean=nanmean(Z2);
xmed = nanmedian(Z2);
xstd = nanstd(Z2);

% String X1 out as a col vector
Z3=Z1';
z3=Z3(:);  % first nfiles rows are seasons 1, 2, ....
% for year 1, next nfiles rows are seasons 1, 2, ... for year 2, and so on

nobs=length(z3); % number of observations (years  times seasons)

% Make a cv , same length as z3, with the seasona means
Y=repmat(xmean,nplot,1);
Y=Y';
ymean=Y(:);

% Subtract monthly means
z3=z3-ymean;


% Compute runs statistics
ksign=-1; % interested in drought, not wetness
xc=0;  % threshold is a negative anomaly
obs=(1:nobs)';
[Px,nx,sx]=runs(z3,obs,xc,ksign);

% Convert sequential seasonal end of runs to year,season
pend=Px(:,2);  % ending sequential obs of drought
yearend = yrplot(1) + floor(pend/nfiles);
seasend = rem(pend,nfiles);
L1temp=seasend==0;
seasend(L1temp)=nfiles;
yearend(L1temp)=yrplot(1)+pend(L1temp)/nfiles-1;


%******************* Figure 1 plot of seasonal info
figure(1)

% Text info summarizing runs with largest run sum
[sum1,isum1]=max(sx);
Ltemp=sx==sum1; % any ties with sum1?
ftemp=find(Ltemp); % indices to multiple maximums
[temp,itemp]=min(nx(ftemp));
isum1=ftemp(itemp);
sum1=sx(isum1);
long1=nx(isum1);  % run-length of run with greatest run-sum
str2a=sprintf('%g, ending in seas#%2.0f, %4.0f',sum1,seasend(isum1),yearend(isum1));
str3a=sprintf(' (L = %4.0f)',long1);

% Text info summarizing runs with largest run length
[num1,inum1]=max(nx);
Ltemp=nx==num1;
ftemp=find(Ltemp);
[temp,itemp]=max(sx(ftemp));
inum1=ftemp(itemp);
num1=nx(inum1);
sum2=sx(inum1); % run-sum of run with longest length
str2b=sprintf('%g seasons, ending in seas#%2.0f, %4.0f',num1,seasend(inum1),yearend(inum1));
str3b=sprintf(' (S = %g %s)',sum2,units);



% Make separate series for those above zero and below zero
z3hi = z3;
z3hi(z3<=0)=NaN;
h1=bar(z3hi);
set(h1,'FaceColor',colors(1,:));
hold on

z3lo =z3;
z3lo(z3>=0)=NaN;
h2=bar(z3lo);
set(h2,'FaceColor',colors(2,:));

% Build y label
ylab1=sprintf(' %4.0f-%4.0f mean (in)',YRS(2,1),YRS(2,2));
ylab2='Departure from ';
ylab3=[ylab2 ylab1];
ylabel(['Departure (' units ')']);
xlabel('Time')
set(gca,'Xtick',[1:nfiles:nobs],'XtickLabel','','XGrid','on');

% Annotate with year
ybot=get(gca,'YLim');
ybot=ybot(1);
for n = 1:nplot;
   text(n*nfiles-(nfiles-1),ybot-.20,int2str(yrplot(n)));
end

% Add info on run sum and run length
pltext(.3,.95,8,['S_{max} = ' str2a str3a]);
pltext(.3,.90,8,['L_{max} = ' str2b str3b]);

% Title
titlab1=sprintf(' %4.0f-%4.0f',YRS(2,1),YRS(2,2));
title(['Seasonal ' dtype ' Departures' ' from ' titlab1 ' Mean'])

hold off

%************************  PLOT OF SEASONAL MEANS

figure(2);
h3=bar(xmean)
title(['Seasonal Means ' dtype ', ' titlab1]);
ylabel([dtype '(' units ')']);
xlabel('Season');
set(gca,'XTickLabel',str2mat(char(cellstr(seasons))));
set(h3,'FaceColor',[1 1 1],'EdgeColor','flat');
