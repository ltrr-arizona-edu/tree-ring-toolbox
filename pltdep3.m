function pltdep3(pfs,jcols,kopt,YRS,colors,seasons)
% pltdep3:  Season-by-season bar charts of time series of departures, with runs summary   
% CALL: pltdep3(pfs,jcols,kopt,YRS,colors,seasons);
%
% Meko 4-28-97
%
%****************** IN ***********************************
%
% pfs{?)}c paths & filenames of .mat files holding seasonalized climate series
%   Example: {'d:\outfls\coolbig.mat','d:\outfls\warmbig.mat'}... for 2-season
% jcols (1 x ?)i;  which column in the files pfs the data should come from
%   Example: [2 2] means get coolbig data from col 2 of first file in pf1, and
%   warmbig from col 2 of second file in pfs
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
% In figure 1,...nseas bar charts of the anomalies from long-term seasonal means
%   for each season, with annotated runs summary
%
%
%*********************** NOTES *************************************
%
% Horizontal line is plotted at zero, indicating departure zero

close('all'); % close figure windows if any are open

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

a=NaN;

% Initialize  matrices to hold the seasonalized data for the nfiles
yrstats=(YRS(2,1):YRS(2,2))'; % year vector for stats period
nstats=length(yrstats);
yrplot=(YRS(1,1):YRS(1,2))';  % for plot period
nplot=length(yrplot);
Z1=a(ones(nplot,1),ones(nfiles,1));
Z2=a(ones(nstats,1),ones(nfiles,1));


%**************  LOOP OVER SEASONAL FILES
clc
disp('Looping over seasons');
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
   
   % Get monthly 'long term' mean and std dev
   xmean=nanmean(z2);
   xstd = nanstd(z2);
   
   % Rename cv of plot series for convenience
   z3=z1;  
   
   nobs=length(z3); % number of years of data in analysis period
   
   % Form departures from mean
   z3=z3-xmean;
   
   % Compute runs statistics
   ksign=-1; % interested in drought, not wetness
   xc=0;  % threshold is a negative anomaly from zero
   [Px,nx,sx]=runs(z3,yrplot,xc,ksign);
   
   % Get ending year of runs
   yearend = Px(:,2);
   
   %******************* Start the figure
   figure(n)
   
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
  % ybot=get(gca,'YLim');
  % ybot=ybot(1);
  % for n = 1:nplot;
 %     text(n*nfiles-(nfiles-1),ybot-.20,int2str(yrplot(n)));
  % end
   
   % Add info on run sum and run length
   pltext(.3,.95,8,['S_{max} = ' str2a str3a]);
   pltext(.3,.90,8,['L_{max} = ' str2b str3b]);
   
   % Title
   titlab1=sprintf(' %4.0f-%4.0f',YRS(2,1),YRS(2,2));
   title([char(seasnm) ' ' dtype ' Departures' ' from ' titlab1 ' Mean']);
   
   hold off
end
