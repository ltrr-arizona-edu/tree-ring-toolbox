function [yearend,monend,nx,sx]=pltdep1a(pf1,jsers,kopt,YRS,thresh,plotinf)
% pltdep1a:  Monthly bar chart of time series of departures, with runs summary
% CALL: pltdep1a(pf1,jsers,kopt,YRS,plotinf);
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
%   kopt(3) threshold option
%     ==1 long term means (and thresh may be empty)
%     ==2 long term medians (ditto)
%     ==3 specified % of long-term mean (thresh is rv)
%     ==4 specified % of long-term median (ditto)
%     ==5 specified departure from long-term mean (ditto)
%     ==6 specified departure from long-term median (ditto)
%     ==7 specified percentile of monthly value (thresh must be scalar)
%     ==8 spefified values as threshold for each month (rv length 12)
%   kopt(4)  option for event being below or above threshold
%     ==1  event defined as value below threshold (e.g., drought using ppt)
%     ==2  event defined as value above threshold
% 
% YRS (1 x 2)i>  start, end year of:
%   row 1: years of data to produce bar chart for -- can be a subset of years of full data
%   row 2: years to compute monthly means and std devs on
% thresh (1 x 12)r (or empty) -- monthly thresholds
%    If kopt(3) is 1 or 2, thresh may be left []
%    If kopt(3) 3-7, thresh is a rv of length 12, each corresponding to a month. For example,
%       [30 30 30 30 30 30 30 30 30 30 30 30] sets threhold at 30th percentile  if
%       kopt(3) is 7. If kopt(3)==8, thresh is rv of specified thresholds
% plotinf{} info for coloring and labeling plots
%    {1} colors(2 x 3)i  colors for upper and lower anomalies
%    {2} ylims (1 x 2)r  y-axis limits (program will allow change in dialog)
%    {3} sernam (1 x ?)s   short (8 chars or so) identifying name for series
%		{4} units (1 x ?)s  string units of data, for labeling y axis (e.g., in)
%		{5} titln2 (1 x ?)s threshold info
%     {6} fmt4 (1 x ?)s  format for printing thresholds in title.
%           Example '%4.1f  '  
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
if length(kopt)~=4;
   error('kopt should be 1 x 4');
end
if kopt(1)==1; % precip
   dtype='PPT';
elseif kopt(1)==2;
   dtype='TMP';
else
   error('kopt(1) must be 1 or 2');
end
if kopt(3)<1 | kopt(3)>8;
   error('kopt(3) must be between 1 and 8, inclusive');
end
if kopt(4)<1 |kopt(4)>2,
   error('Invalid kopt(4)')
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
units=plotinf{4}; % units of data
titln2=plotinf{5}; % line 2 of title -- threshold logic
fmt4=plotinf{6}; % format for a threshold value in title

% thresh
if kopt(3)>2 & kopt(3)~=7;
   if length(thresh)~=12;
      error('thresh must be rv of length 12 for this kopt(3)');
   end
elseif kopt(3)==1 | kopt(3)==2; % 
   if ~isempty(thresh) & length(thresh)~=12;
      error('thresh must be [] or 12-element rv for this kopt(3)');
   end
else; % kopt(3)==7
   if length(thresh)~=1;
      error('thresh must be scalar percentile for this kopt(3)');
   end
end


% Make pointer to rows of X
L2 = yr>=YRS(2,1) & yr<=YRS(2,2); % years to compute means and std devs on
L1 = yr>=YRS(1,1) & yr<=YRS(1,2); % years to plot

% Get monthly 'long term' meansm medians and std devs
xmean=nanmean(X(L2,:));
xmed=nanmedian(X(L2,:));
xstd = nanstd(X(L2,:));

% Get the subset of rows (years)of X to be plotted
X1 = X(L1,:);
[m1,n1]=size(X1);
yr2=yr(L1,:);

% String X1 out as a col vector
X2=X1';
x2=X2(:);  % first 12 rows are jan-dec year 1, next 12 rows are jan-dec year 2, etc
nmonths=length(x2);

% Make a cv , same length as x2, with the monthly thresholds
if kopt(3)==1; % threshold is monthly mean
   Y=repmat(xmean,m1,1);
   xtemp=xmean;
elseif kopt(3)==2; % medians
   Y=repmat(xmed,m1,1);
   xtemp=xmed;
elseif kopt(3)==3; % percent of long-term mean
   xtemp=xmean .* thresh/100;
   Y=repmat(xtemp,m1,1);
elseif kopt(3)==4;  % percent of long term median
   xtemp=xmed .* thresh/100;
   Y=repmat(xtemp,m1,1);
elseif kopt(3)==5; % threshold is mean plus specified amt 
   xtemp=xmean + thresh;
   Y=repmat(xtemp,m1,1);
elseif kopt(3)==5;  % threshold is median plus some amount
   xtemp=xmed + thresh;
   Y=repmat(xtemp,m1,1);
else; % kopt(3)==7, and using percentile
   xtemp=prctile(X(L2,:),thresh);
   Y=repmat(xtemp,m1,1);
end
   
 
Y=Y';
ybase=Y(:); % baseline from which departures are calculated

% Subtract monthly baseline
x2=x2-ybase;


% Compute runs statistics
if kopt(4)==1;  % event is value below threshold
   ksign=-1;
else
   ksign=1;
end

xc=0;  % baseline of converted series for for computing departures is zero
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
   if rem(yr2(n),2)==0;
      hbot(n)=text(n*12-12,ybot,int2str(yr2(n)),...
         'HorizontalAlignment','center','VerticalAlignment','top');

   end
   
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
   % Annotate with year
   ybot=get(gca,'YLim');
   ybot=ybot(1); 
   for n = 1:m1;
      if rem(yr2(n),2)==0;
         oldpos=get(hbot(n),'Position');
         oldpos(2)=ybot;
         set(hbot(n),'Position',oldpos);
      end
   end
   
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
titln1=[sernam ':  Monthly ' dtype ' Departures, '  titlab1 ];
titln3=sprintf(fmt4,xtemp(1:length(xtemp)));
title({titln1,titln2,titln3});

set(gca,'Position',[0.1300    0.1100    0.7750    0.77]);

hold off

%************************  PLOT OF MONTHLY MEANS WITH THRESHOLDS

figure(2);
barmeat=[xmean' xtemp'];
bar(barmeat)
title('Long Term Monthly Means and Thresholds');
ylabel([dtype '(' units ')']);
xlabel('Month');
legend('Means','Thresholds');
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