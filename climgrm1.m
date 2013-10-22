function climgrm1(P,T,txtin,windno,colors,axlim)
% climgrm1: a single McClaran climagram
% CALL:  climgrm1(P,T,txtin,windno,colors,axlim);
%
% Meko 7-7-99
%
%***************** IN 
%
% P (1 x 12)r  monthly means of precipitation
% T (1 x 12)r  monthly means of temperature
% txtin{} -- Labeling for plots
%   {1} (1 x ?)s title for plot (will go above top plot). Say,
%       'Climagram for Rapid City, Period 1901-99'
%   {2} (1 x ?)s y axis label for pcp  Example
%      'Precipitation (mm)'
%   {3} (1 x ?)s y axis label for tmp
% windno: window number for plot
% colors (2 x 3)r  colors for bars, or ppt (row 1) and lines, or temp (row 2)
%     [] means accept default
% axlim (1 x 4)r limits for left & right y axes
%      [LowerL UpperL  LowerR UpperR]
%     [] means let data dictate
%
%
%********  NOTE
%
% Makes two axes in a figure window.  Top axes is boxplot distribution of monthly
% precip.  Bottom axes is for monthly temperature, which might be monthly means of
% daily means, daily maximums, or whatever.  
%
% Labeling kept to a minimum to allow flexibility of user to label with correct
% units in calling program
%
% Be sure to have no existing figure window windno  upon calling climgrm1

%------ CHECK INPUT
[mtemp,ntemp]=size(P);
if ntemp~=12 |mtemp~=1;
   error('P must be 1 x12');
end
[mtemp,ntemp]=size(T);
if ntemp~=12| mtemp ~=1;
   error('T must be 1 x 12');
end


%************ FIGURE
% Caption
% Figure .  Graphs showing distributions of monthly climatic variables at
% Great Sand Dunes.  Top: monthly total precipitation.  Bottom: monthly averaged
% daily maximum temperature.  Boxplots based on data from Great Sand Dunes National 
% Monument. Each box summarizes distribution of 47 values (1951-97). Boxplot 
% elements are median (horizontal line in middle of box); interquartile range 
% (extents of box); adjacent values (whiskers), defined as the most extreme values
% still within 1.5 times the interquartile range of the edges of the box; and 
% outliers (+) -- values greater than that distance from the edges of the box. 
% Boxplot definitions from Cleaveland (1993). 
%close all;
%figure(windno);


set(gcf,'DefaultLineLineWidth',1.5);
t = (1:12);
xtlab = {'J','F','M','A','M','J','J','A','S','O','N','D'};

fun1='bar'; fun2='plot';
[ax,h1,h2]=plotyy(t,P,t,T,fun1,fun2);

% If specify y-axis limits, will also want uniform y-axis ticking and labeling
if ~isempty(axlim);
   PTick=axlim(1):1:axlim(2); % sets at each inch
   TTick=axlim(3):20:axlim(4); % every 20 deg
end;


% There is a Matlab bug that results in ticks being placed on both y axes when
% you try to set the YTick property of the left y-axis. The following is
% a workaround
hlab1=get(ax(1),'YLabel');
set(hlab1,'string',txtin{2});

% Mitch wants horizontal yaxis label on right
set(hlab1,'Rotation',0);
set(hlab1,'Position',[-1.5000    2.4780   17.3205]);
set(hlab1,'HorizontalAlignment','Right');

set(ax(1),'XTick',[1:12],'XTickLabel',xtlab);
if ~isempty(axlim);
   set(ax(1),'YLim',axlim(1:2));
   set(ax(1),'YTick',PTick);
   set(ax(1),'XLim',[-1 13]);
end;

hlab2=get(ax(2),'YLabel');
set(hlab2,'string',txtin{3});

% Mitch wants horizontal y-axis labels on right, too
set(hlab2,'Rotation',0);
set(hlab2,'HorizontalAlignment','Right');
set(hlab2,'Position',[14.5000   70.0000   17.3205]);

set(ax(2),'XTick',[]);
if ~isempty(axlim);
   set(ax(2),'YLim',axlim(3:4));
   set(ax(2),'YTick',TTick);
   set(ax(2),'XLim',[-1 13]);
     
end;
set(h2,'LineWidth',2);

if ~isempty(colors);
   set(h1,'FaceColor',colors(1,:));
   set(h1,'EdgeColor',colors(1,:));
   set(h2,'Color',colors(2,:));
   set(ax(1),'Ycolor',[0 0 0]);
   set(ax(2),'Ycolor',[0 0 0]);
end

% Expand font size
set(ax(1),'FontSize',12,'box','off');
set(ax(2),'FontSize',12,'box','off');

title(txtin{1},'FontSize',12);
disp('here')
