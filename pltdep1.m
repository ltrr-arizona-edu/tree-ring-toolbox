function [yearend,monend,nx,sx]=pltdep1(pf1,jsers,kopt,YRS,plotinf)
% pltdep1:  Monthly bar chart of time series of departures, with runs summary
% CALL: pltdep1(pf1,jsers,kopt,YRS,plotinf);
%
% Meko 9-28-97
%
%****************** IN ***********************************
%
% pf1(1 x ?)c path & filename of .mat file holding monthly clim data
%   13-col year plus monthly data assumed to be in matrix named X or Y
% jsers (1 x 1)i  seq number of desired series in X or Y (see notes)
% kopt(1 x ?)i  options
%   kopt(1) ==1 ppt
%           ==2 temperature
%   kopt(2) ==1 plot anomalies in original units
%           ==2 plot standardized anomalies (not implemented yet)
% YRS (1 x 2)i>  start, end year of:
%   row 1: years of data to produce bar chart for -- can be a subset of years of full data
%   row 2: years to compute monthly means and std devs on
% plotinf{} info for coloring and labeling plots
%    {1} colors(2 x 3)i  colors for upper and lower anomalies
%    {2} ylims (1 x 2)r  y-axis limits (program will allow change in dialog)
%    {3} sernam (1 x ?)s   short (8 chars or so) identifying name for series
%
%******************** OUT **************************************
%
% In figure 1, a bar chart of the anomalies from long-term monthly means
% In figure 2, a bar chart of long-term monthly means
%
%
%*********************** NOTES *************************************
% jsers -- the input matrix (X or Y) is assumed to contain a year column and
%   groups of 13 columns of monthly climate data for one or more stations.  jsers
%   tells which of the stations (or regions, gridpoints, or whatever) you want to
%   use.  For example, say rega.mat has Z, with 25 cols representing monthly data
%   for N san pedro and S san pedro.  jsers==1 says pull off cols 2:13, which is
%   data for N san pedro
% Horizontal line is plotted at zero, indicating departure zero


% Get the .mat file that has matrix
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

% Check that X has correct number of cols
[mX,nX]=size(X);
yr = X(:,1);
if rem(nX,12)~=1;  
   error('X does not have year col plus multiple of 12 data cols')
end

% Lop year col off X
X(:,1)=[];

% Compute the col indices into X for desired series, then cull the desired cols
jtemp=[1:12]+12*(jsers-1);
X=X(:,jtemp);

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


% Check out YRS
yrgo1=yr(1);
yrsp1=yr(mX);
[mtemp,ntemp]=size(YRS);
if ~(mtemp==2 & ntemp==2);
   error('YRS should be 2 x 2');
end
if any(YRS(:,1)<yrgo1);
   error('Starting years in YRS inconsistent with years in data matrix');
end
if any(YRS(:,2)>yrsp1);
   error('End years in YRS inconsistent with years in data matrix');
end

%-------  plotinf
colors=plotinf{1}; % colors for bars
ylims=plotinf{2}; % limits for y axis
sernam=plotinf{3}; % Series name -- used in plot title


% Make pointer to rows of X
L2 = yr>=YRS(2,1) & yr<=YRS(2,2); % years to compute means and std devs on
L1 = yr>=YRS(1,1) & yr<=YRS(1,2); % years to plot

% Get monthly 'long term' means and std devs
xmean=nanmean(X(L2,:));
xstd = nanstd(X(L2,:));

% Get the subset of rows (years)of X to be plotted
X1 = X(L1,:);
[m1,n1]=size(X1);
yr2=yr(L1,:);

% String X1 out as a col vector
X2=X1';
x2=X2(:);  % first 12 rows are jan-dec year 1, next 12 rows are jan-dec year 2, etc
nmonths=length(x2);

% Make a cv , same length as x2, with the monthly means
Y=repmat(xmean,m1,1);
Y=Y';
ymean=Y(:);

% Subtract monthly means
x2=x2-ymean;


% Compute runs statistics
ksign=-1; % interested in drought, not wetness
xc=0;  % threshold is a negative anomaly
mons=(1:nmonths)';
[Px,nx,sx]=runs(x2,mons,xc,ksign);

% Convert sequential month of beginning of runs to year,month
pend=Px(:,2);  % ending sequential obs of drought
monend=rem(pend,12);
L1temp=monend==0;
yearend=yr2(1)+floor(pend/12);
monend(L1temp)=12;
yearend(L1temp)=yr2(1)+pend(L1temp)/12-1;

%if monend==0;
  % monend=12
   %yearend=yr2(1)+pend/12-1;
%else
   %yearend = yr2(1) + floor(pend/12);
   %monend = rem(pend,12);
%end


% Text info summarizing runs with largest run sum
[sum1,isum1]=max(sx);
Ltemp=sx==sum1; % any ties with sum1?
ftemp=find(Ltemp); % indices to multiple maximums
[temp,itemp]=min(nx(ftemp));
isum1=ftemp(itemp);
sum1=sx(isum1);
long1=nx(isum1);  % run-length of run with greatest run-sum
str2a=sprintf('%g, ending in %2.0f/%4.0f',sum1,monend(isum1),yearend(isum1));
str3a=sprintf(' (L = %4.0f)',long1);


% Text info summarizing runs with largest run length
[num1,inum1]=max(nx);
Ltemp=nx==num1;
ftemp=find(Ltemp);
[temp,itemp]=max(sx(ftemp));
inum1=ftemp(itemp);
num1=nx(inum1);
sum2=sx(inum1); % run-sum of run with longest length
str2b=sprintf('%g months, ending in %2.0f/%4.0f',num1,monend(inum1),yearend(inum1));
str3b=sprintf(' (S = %g %s)',sum2,units);


% Make separate series for those above zero and below zero
x2hi = x2;
x2hi(x2<=0)=NaN;
h1=bar(x2hi);
set(h1,'FaceColor',colors(1,:));
hold on

x2lo =x2;
x2lo(x2>=0)=NaN;
h2=bar(x2lo);
set(h2,'FaceColor',colors(2,:));

% Build y label
ylabel(['Departure (' units ')']);
% xlabel('Time')
set(gca,'Xtick',[1:12:nmonths],'XtickLabel','','XGrid','on');

% Annotate with year
ybot=get(gca,'YLim');
ybot=ybot(1);
for n = 1:m1;
   text(n*12-9,ybot-.15,int2str(yr2(n)));
end


% Optionally change limits of y axis

ktemp1=questdlg('Want control over y-axis limits?');
switch ktemp1;
case 'No'
case 'Cancel';
case 'Yes';
   ktemp2=questdlg('Use the hard-coded limits passed in plotinf{2}?');
   switch ktemp2;
   case {'No','Cancel'};
      prompt1={'Enter max and min value for y axis'};
      title1='OPTIONAL CONTROL OF Y-AXIS LIMITS';
      lineno=1;
      defans=get(gca,'YLim');
      ans1=dlgi001(title1,prompt1,lineno,defans)
   case 'Yes';
      ans1=ylims;
   end
   set(gca,'YLim',ans1);
end




% Optionally, add info on run sum and run length
ktemp=questdlg('Want to annotate plot with highest run sum and run length?')
switch ktemp
case 'No';
case 'Cancel';
case 'Yes'
   pltext(.3,.95,8,['S_{max} = ' str2a str3a]);
   pltext(.3,.90,8,['L_{max} = ' str2b str3b]);
end


% Title
titlab1=sprintf('%4.0f-%4.0f',YRS(1,1),YRS(1,2));
titlab2=sprintf(' %4.0f-%4.0f',YRS(2,1),YRS(2,2));
title([sernam ':  Monthly ' dtype ' Departures, '  titlab1 ',  from ' titlab2 ' Mean'])

hold off

%************************  PLOT OF MONTHLY MEANS

figure(2);
bar(xmean)
title(['Monthly Mean ' dtype ', ' titlab1]);
ylabel([dtype '(' units ')']);
xlabel('Month');
set(gca,'XTickLabel',['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'])

%************************* ARRANGE Px,nx,sx in decr order of run sum

[ssx,ssxi]=sort(sx);
sx=flipud(sx(ssxi));
nx=flipud(nx(ssxi));
yearend=flipud(yearend(ssxi,:));
monend=flipud(monend(ssxi,:));

%**************************** MAKE ASCII FILE OF RESULTS, WITH
% RUNS ORDERED IN DESCENDING ORDER OF RUN-SUM
ktemp=questdlg('Want an ascii file summarizing runs?')
switch ktemp
case 'No';
case 'Cancel';
case 'Yes'
   [file2,path2]=uiputfile('*.txt','Ascii summary table of runs');
   pf2=[path2 file2];
   fid2=fopen(pf2,'w');
   head1='Runs Ranked in Descending order of Run Sum';
   head2=['Series = ' sernam];
   head3=['Analysis period: ' titlab1];
   head4=['   Long Term Means from: ' titlab2];
   
   head5=' Rank  End Month  Sum   Length(mos) '
   nlist=(1:length(nx))';
   
   % Allow for making ascii table not have the full list of runs
   prompt1={'Enter number of runs to list'};
   title1='';
   lineno=1;
   defans={size(sx,1)};
   ans1=inputdlg(prompt1,title1,lineno,defans);
   ans1=str2num(ans1{1});
   
   Z=[nlist monend yearend sx nx];
   Z=Z(1:ans1,:);
   
   fmt1='%3.0f \t%2.0f/%4.0f \t%g\t(%2.0f)\n';
   fprintf(fid2,'%s\n%s\n%s.   %s\n\n',head1,head2,head3,head4);
   fprintf(fid2,'%s\n\n',head5);
   fprintf(fid2,fmt1,Z');
   fclose(fid2);
end

