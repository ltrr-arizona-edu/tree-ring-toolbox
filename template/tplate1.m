% tplate1.m  script graphic with automatic prep for web
%
% Background color, and properties for lines, text, arrows are changed depending on medium of output.
% Web effects call for black background, correspondingly light colors, larger text, etc.
% Colors are not hard coded.  Use natural color order instead.


% Clear variables and close any windows
clear;
close all;

% Hard code data directory
dir1 = 'c:\projs\ai3\sacrflow\';
pf1=[dir1 'rfog3'];


% Prompt for background of figures
kmen1=menu('Choose color for figure background',...
   'Black -- best for on-screen viewing',...
   'White -- best for printout');
if kmen1==1;
   bkgcolor='Black';
else;
   bkgcolor='White';
end;
switch bkgcolor;
case 'Black';
   whitebg([0 0 0]);
   colnull=[1 1 1];
case 'White';
   whitebg([1 1 1]);
   colnull=[0 0 0];
end;
close all;

% Get data and generate plot
eval(['load ' pf1 ' Y2 y yry;']);
ymean=mean(exp(log(10)*y));
L=~isnan(Y2(:,2));
Y2=Y2(L,:);
yrx=[Y2(:,1) ; yry];
x=[Y2(:,2) ;  y];
clear y yry;
x=exp(log(10)*x);

figure(1);
b=wtsgaus(100);
[y,yry]=mafilt1(x,yrx,94,2);
%[y,yry]=filter1(x,yrx,b,1);
L=yry==1112 | yry==1350;
i=find(L);
hp1=plot(yry,y,yry(i),y(i),'o',[min(yry) max(yry)],[ymean ymean]);
xlabel('Ending Year of 94-Year Period');
ylabel('Flow (MAF)');
title('94-Year Running Mean of Reconstructed Flow');
cord=get(gca,'ColorOrder');
set(gca,'YLim',[15 23]);

% Arrows & text
text(1000*1.2369,1000*0.0197,'A.D. 1112','HorizontalAlignment','Center');
text(1000*1.4802,1000*0.0197,'A.D. 1350','HorizontalAlignment','Center');
harr1=arrow(1000*[1.2147    0.0194],1000*[1.1290    0.0183]);
harr2=arrow(1000*[1.4636    0.0194],1000*[1.3641    0.0176]);

% TAILOR OUTPUT FEATURES TO MEDIUM
switch bkgcolor;
case 'Black';
   fsbase=get(gca,'FontSize');
   
   % Font size for title, axes, and axes labels
   set(gca,'FontSize',fsbase+2);
   % Title font
   htemp=get(gca,'Title');
   set(htemp,'FontSize',fsbase+4);
   % Axes label font
   htemp=get(gca,'XLabel');
   set(htemp,'FontSize',fsbase+2);
   htemp=get(gca,'YLabel');
   set(htemp,'FontSize',fsbase+2);
   
   % Line widths for lines
   htemp=findobj('Type','Line');
   nline=length(htemp);
   for j=1:nline;
      lwthis=get(htemp(j),'LineWidth');
      set(htemp(j),'LineWidth',lwthis+1);
   end;
   
   % Text
   htext=findobj('Tag','Text');
   if ~isempty(htext);
      ntext=length(htext);
      for j=1:ntext;
         fstemp=get(htext(j),'FontSize');
         set(htext(j),'FontSize',fstemp+2);
      end;
   end;
      
   % Arrows: line width and face and edge colors
   harrow=findobj('Tag','Arrow');
   if ~isempty(harrow);
      narrow=length(harrow);
      for j=1:narrow;
         lwtemp=get(harrow(j),'LineWidth');
         set(harrow(j),'Edgecolor',[1 1 1]);
         set(harrow(j),'Facecolor',[1 1 1]);
         set(harrow(j),'LineWidth',lwtemp+1);

      end;
   end;
   
        
case 'White';
otherwise;
end;


grid;
zoom xon;

figure(2);
f=[0:.001:.04]';
u=freqres2(b,1,f);
plot(f,u);
grid;
zoom xon;