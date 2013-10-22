function colmap02(txtin,datin)
%colmap01:  color map tree-ring seasonal response to climate. See resp02.m
%CALL: colmap01(txtin,datin);
%
% Meko 10-28-97
%
%**************  IN **************************************
%
% txtin{} -- text info for figure
%   {1} tit(1 x ?)s  title
%   {2} xlab(1 x ?)s  xlabel
%   {3} ylab(1 x ?)s  ylabel 
%   {4} cmap1(4 x 3)r  user designed color map to use
% datin{} -- data
%   {1} endmo (1 x 1)i end month of the 'tree' year (e.g., 8 == august)
%   {2} Z (12 x 12)r  matrix of data to be plotted
%   {3} kopt(1 x 1)i  options
%      kopt(1):  indexed color or color values
%         ==1 color values computed from data D
%         ==2 indexes to discrete set of color computed from D (not
%					implemented yet)
%   {4} clim (1 x 2)r  color limits for mapping
%
%**********  OUT ******************************
%
% No output.  Function just gives a color figure
%
%******************** NOTES **************************
%
% Color mapping.  Optionally, colmap01 makes a new color map with a constant color
% for any specified range of the variable in Z. This by kopt(1).  Otherwise
% colmap01 interpolates within the color map for a continuum of colors

a=NaN;

%***************  UNLOAD CELLS

endmo=datin{1};
Z=datin{2};
kopt(1)=datin{3};
clim=datin{4};

tit=txtin{1};
xlab=txtin{2};
ylab=txtin{3};
cmap1=txtin{4};





%*************** Add null 13th row and column to D because  pcolor
% does not use them
Z1=[Z a(ones(12,1),:)];
Z1=[Z1; a(:,ones(13,1))];


%**************** Color map

colormap(cmap1);


pcolor(Z1);


% --------------- axes labels

xlabel(xlab);
ylabel(ylab);

set(gca,'Xtick',[2.5 4.5 6.5 8.5 10.5 12.5],...
   'XtickLabel',[' 2';' 4';' 6';' 8';'10';'12']);
set(gca,'Ytick',[1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 12.5]);

set(gca,'Clim',clim);

% Must build y tick labels from month letters, with correct ending month
mons24 = 'JFMAMJJASONDJFMAMJJASOND'; % 24 months
j2 = 12+endmo;
j1 = j2-11;
mons12 = (mons24(j1:j2))';
set(gca,'YTickLabel',mons12);

title(tit)


ll=[13.5 4.5; 13.5 7.5];
ur=[14.0 5.0; 14.0 8.0];
txtlab=['+';'-'];
fntsz=10;      

h2=colorbar;
set(h2,'Ytick',[11 55],...
   'YtickLabel',['-';'+'],...
   'Fontsize',24)


