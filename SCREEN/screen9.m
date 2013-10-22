function screen9
% screen9a: map of site locations -- original, time screened, ppt screened
% CALL: screen9a
%
% Meko 8-19-97
%
%******************  IN **********************************
%
% These files must be accessible:
%
%   ap1os5.dat -- ascii all-numeric file with information on the final modeling
%     This file has 20 cols, only some of which are needed
%     col 1,2 long and lat in dec degree plotting format
%   gospa.dat -- start, end year of chronolgies.  Needed to rebuild pointer for
%     year-coverage screening
%   ap1os6.mat -- holds logical pointers Lcull and Lneg to full-set chrons
%     telling which have been selected (Lcull) and which have negative r with
%     the climate variable. Also holds yrgo and yrsp, the start and end years
%     by which chrons were initially screened befor climate screening

% Close any maps or figure windows
close all



% Get the data output by screen5.m
[file1,path1]=uigetfile('*os5.dat','Infile of data output by screen5.m');
pf1=[path1 file1];
eval(['load ' pf1]);
eval (['Z = ' lower(strtok(file1,'.')) ';']);


% Pull longitudes and latitudes for full set
lonall=Z(:,1); % longitudes
latall=Z(:,2); % latitudes

% Get the start year and stop year info
[file2,path2]=uigetfile('gosp*.dat','Infile of full set start and end yrs');
pf2=[path2 file2];
eval(['load ' pf2]);
eval (['xyall = ' lower(strtok(file2,'.')) ';']);

% Get the logical pointer to final selected chrons
[file3,path3]=uigetfile('*os6.mat','Infile of data output by screen6.m');
pf3=[path3 file3];
eval(['load ' pf3]);  % Puts Lcull, Lneg, logscrn in workspace


% Make a pointer to acceptable year coverage
Lyr = gospa(:,1)<=yrgo & gospa(:,2)>=yrsp;


% Pointer to chrons with negative tree-climate simple ar
% Lneg is this.

% Pointer to final selected chrons
% Lcull is this

% Pointer to OE(1,0) models
%L5=ap1os5(:,7)==1 & ap1os5(:,8)==0;



%***************  Load vector files of longitude and latitude
load worldlo
load usalo


%-----------Make mutually exclusive pointers to the types of points

% black -- full set chrons not with desired year coverage
Lk = ~Lyr;
% Magenta -- sites with OK time coverage, but screened out for ppt signal weak
Lm = Lyr & ~Lcull;
% green -- final selected
Lg = Lcull;


%*****************  MAP 1 -- FULL SET AND SELECTED PPT-SENSITIVE
figure(1);

% Note: full limits 22-70N, -150 to -96W
axesm('MapProjection','eqaconic','Frame','on',...
	'MapLatLimit',[22 70],...
	'MapLonLimit',[-150 -96]);
gridm('MLineLocation',10,'PLineLocation',10,'GLineStyle','-');
mlabel('MLabelLocation',10,'MLabelParallel','south')
plabel('PLabelLocation',10);
displaym(POline);
plotm(statelat,statelon);

% Plot data points

h1=plotm(latall(Lk),lonall(Lk),'m.')
h2=plotm(latall(Lm),lonall(Lm),'k.');

%plotm(latall(Lm),lonall(Lm),'m.');
set(gcf,'Position',[0 120 350 420]);
set(gca,'position',[0 0.05 1 .90]);

set(gcf,'PaperOrientation','Portrait');
%set(gcf,'PaperPosition',[.25 .25 8.0 10.5]);  % for 'orient tall' default
set(gcf,'PaperPosition',[.25 .25 4.0 5]);  
txt1 = sprintf('%3.0f of %3.0f chrons covering %4.0f-%4.0f',...
      sum(Lyr),length(Lcull),yrgo,yrsp);
title(txt1);
%h=[h1 h2];
%legend(h,'Full 669 site network','Subset covering 1700-1979');



%*****************  MAP 2 -- year screened  AND SELECTED PPT-SENSITIVE
figure(2);

% Note: full limits 22-70N, -150 to -96W
axesm('MapProjection','eqaconic','Frame','on',...
	'MapLatLimit',[22 70],...
	'MapLonLimit',[-150 -96]);
gridm('MLineLocation',10,'PLineLocation',10,'GLineStyle','-');
mlabel('MLabelLocation',10,'MLabelParallel','south')
plabel('PLabelLocation',10);
displaym(POline);
plotm(statelat,statelon);

% Plot data points

h1=plotm(latall(Lm),lonall(Lm),'g.');
h2=plotm(latall(Lg),lonall(Lg),'k.');
%h=[h1 h2];

%plotm(latall(Lm),lonall(Lm),'m.');
set(gcf,'Position',[0 120 350 420]);
set(gca,'position',[0 0.05 1 .90]);

set(gcf,'PaperOrientation','Portrait');
%set(gcf,'PaperPosition',[.25 .25 8.0 10.5]);  % for 'orient tall' default
set(gcf,'PaperPosition',[.25 .25 4.0 5]);  
%legend(h,'88/343 chrons ppt-sensitive');
txt2=sprintf('%3.0f of %3.0f chrons selected',sum(Lcull),sum(Lyr))
title(txt2);

