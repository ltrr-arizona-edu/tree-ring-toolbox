function [mpcp,npcp,mtmp,ntmp]=climgram(P,T,txtin,windno)
% climgram: graph monthly distribution of precipitation and temperature
% CALL:  [mpcp,npcp,mtmp,ntmp]=climgram(P,T,txtin,windno);
%
% Meko 7-2-98
%
%***************** IN 
%
% P (? x 13)r  monthly precipitation matrix, year in col 1
% T (? x 13)r  monthly temperature matrix, year in col 1
% txtin{} -- Labeling for plots
%   {1} (1 x ?)s title for plot (will go above top plot). Say,
%       'Climagrams for Rapid City'
%   {2} (1 x ?)s y axis label for pcp  Example
%      'Precipitation (mm)'
%   {3} (1 x ?)s y axis label for tmp
% windno: window number for plot
%
%***************** OUT
%
% mpcp (1 x 12)r  means for pcp
% npcp (1 x 12)n  number of years in sample
% mtmp (1 x 12)r  means for tmp
% ntmp (1 x 12)n  number of years in sample
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
% Be sure to have no existing figure window windno  upon calling climgram.m

%------ CHECK INPUT
[mtemp,ntemp]=size(P);
if ntemp~=13;
   error('P must have 13 cols');
end
[mtemp,ntemp]=size(T);
if ntemp~=13;
   error('T must have 13 cols');
end

%------ Compute mean of monthly precip, using only valid data
P1 = P(:,2:13);
npcp = sum(~isnan(P1));
mpcp = nanmean(P1);

%------ Compute mean of monthly temperature, using only valid data
T1 = T(:,2:13);
ntmp = sum(~isnan(T1));
mtmp = nanmean(T1);

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

figure(windno);
set(gcf,'DefaultLineLineWidth',1.5);
t = (1:12);
xtlab = {'J','F','M','A','M','J','J','A','S','O','N','D'};

hax1=axes('Position',[.1 .1 .8 .4]);
hax2=axes('Position',[.1 .55 .8 .4]);

%************* TMP
axes(hax1);
boxplot(T1);
set(gca,'XTickLabel',[]);
set(gca,'Xtick',[1:12],'XTickLabel',xtlab);
set(gca,'Fontsize',11);

ylabel([]);
ylabel(txtin{3},'FontSize',14);
xlabel(' ');

%************* PCP
axes(hax2);
boxplot(P1);
set(gca,'XTickLabel',[]);
set(gca,'Xtick',[1:12],'XTickLabel',' ');
ylabel([]);
set(gca,'Fontsize',11);

ylabel(txtin{2},'FontSize',14);
xlabel(' ');
title(txtin{1});
