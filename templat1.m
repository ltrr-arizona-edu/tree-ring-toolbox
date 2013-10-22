% templat1.m
% template for slide vs web vs paper vs poster

% Some iniatial default settings
set(0,'DefaultTextFontName','arial');  % arial font 
MT1='x'; MT2='o'; MT3='*';  MT4='square'; MT5='+'; % markers
MS1=5; MS2=5; MS3=5; MS4=5; MS=5; % marker sizes

% This will ensure that initial menus are readable
colordef white;
set(gcf,'Color',[.8 .8 .8]);

% Figure background colors, potential  
FigC1=[.8 .8 .8]; FigC2=[.2 .2 .2];  FigC3=[0 0 0]; FigC4=[.4 .4 .4];

% Axis background colors, potential. 
AC1=[1 1 1];  AC2=[0 0 0];


%*** CHOOSE TYPE OF MEDIA

kmen1=menu('Choose media',...
   'Paper',...
   'Transparency',...
   'Slide',...
   'Web');
if kmen1==1; % paper
   media='Paper';
elseif kmen1; % transparency
   media='Transp';
elseif kmen1==3; % slide
   media='Slide';
elseif kmen1==4; % web 
   media='Web';
   colordef black;
else;
   error('Invalid media');
end;


%*** INITIAL SETTINGS OF COLOR DEFINITION AND FIGURE COLOR
figure(1);
if kmen1==1; % paper
   FigC=FigC1;
   media='Paper';
   colordef white;
   set(gcf,'Color',FigC);
elseif kmen1==2; % transparency
   FigC=FigC1;
   media='Transp';
   colordef white;
   set(gcf,'Color',FigC);
elseif kmen1==3; % slide
   media='Slide';
   FigC=FigC2;
   colordef black;
   set(gcf,'Color',FigC);
elseif kmen1==4; % web 
   media='Web';
   FigC=FigC2;
   colordef black;
   set(gcf,'Color',FigC);
else;
   error('Invalid media');
end;


%*** INITIAL SETTINGS OF SOME OTHER FIGURE PROPERTIES, AND OF AXIS, AND LINE PROPERTIES
switch media;
case 'Paper' ;
   FSaxis=10;
      kmen2=menu('Choose one','Color','BW');
   if kmen2==1; % color
      AC=[1 1 1]; % axis background
      COL0=[0 0 0];  % axis lines
      COL1=[0 0 1];  % line 1
      COL2=[.5 0 .5];  % line2
      COL3=[0 .5 .5]; COL4=[.5 0 0];  COL5=[0 0 0];
      LT1='-'; LT2='-'; LT3='-'; LT4='-'; LT5='none';
      LW1=2; LW2=2; LW3=2; LW4=2; LW5=2;
   elseif kmen2==2; % bw
      AC=[1 1 1]; % axis background
      COL0=[0 0 0];  % axis lines
      COL2=[0 0 0];  % line 1
      COL1=[0 0 0];  % line2
      COL3=[0 0 0]; COL4=[0 0 0];  COL5=[0 0 0];
      LT1='-'; LT2=':'; LT3='--'; LT4='-.'; LT5='none';
      LW1=0.5; LW2=0.5; LW3=0.5; LW4=0.5; LW5=0.5;
      kmen3=menu('Choose one','Thick gray line2','Thin Black dotted line2');
      if kmen3==1;  % thick gray line 2
         COL2=[.6 .6 .6]; LW2=2;  LT2='-'; 
      end;
      kmen4=menu('Orientation',...
         'Portrait',...
         'Landscape',...
         'Tall');
      if kmen4==1;
         orient='Portrait';
      elseif kmen4==2;
         orient='Landscape';
      elseif kmen4==3;
         orient='Tall';
      else;
         error('Invalid orientation');
      end;
      kmen5=menu('Choose One',...
         'Figure full width of page',...
         'Figure half width of page');
      if kmen5==1;
         figlook='Full-width';
      elseif kmen5==2;
         figlook='Half-width';
         if strcmp(orient,'Landscape');
            error('Half width figlook only makes sense with portrait');
         end;
         
      end;
      
   end;
case 'Transp';
   FSaxis=10;
   AC=[1 1 1]; % axis background
   COL0=[0 0 0];  % axis lines
   COL1=[0 0 1];  % line 1
   COL2=[.5 0 .5];  % line2
   COL3=[0 .5 .5]; COL4=[.5 0 0];  COL5=[0 0 0];
   LT1='-'; LT2='-'; LT3='-'; LT4='-'; LT5='none';
   LW1=2; LW2=2; LW3=2; LW4=2; LW5=2;

case 'Slide';
   FSaxis=14;
   FSaxis=14;
   COL0=[1 1 1];  % axis lines
   COL1=[0 1 0];  % line 1
   COL2=[0 .4 1];  % line2
   COL3=[1 0 0]; COL4=[1 1 0];  COL5=[1 0 1];
   LT1='-'; LT2='-'; LT3='-';
   LW1=2.0;  LW2=2; LW3=2; LW4=2; LW5=2;
   MT1='.'; MT2='o'; MT3='square';  MT4='diamond'; MT5='*';
   MS1=16;   
case 'Web';
   FSaxis=14;
   COL0=[1 1 1];  % axis lines
   COL1=[0 1 0];  % line 1
   COL2=[0 .4 1];  % line2
   COL3=[1 0 0]; COL4=[1 1 0];  COL5=[1 0 1];
   LT1='-'; LT2='-'; LT3='-';
   LW1=2.0;  LW2=2; LW3=2; LW4=2; LW5=2;
   MT1='.'; MT2='o'; MT3='square';  MT4='diamond'; MT5='*';
   MS1=16;   MS2=5; MS3=5; MS4=5; MS5=5;
otherwise;
end;


% SET FONT SIZES FOR AXIS LABELS, TITLE, AND POINT LABELS

FSlab=FSaxis+1;
FStit=FSaxis+3;
FSdot=FSaxis-1;


%********** DUMMY DATA SECTION ****************

load lwwstd.dat;
x1=lwwstd;
load lwwres.dat;
x2=lwwres;
clear lwwstd lwwres;

% Reverse order of plotting series 1 and 2 if thick gray line option
hp1=plot(x1(:,1),x1(:,2),x2(:,1),x2(:,2));
set(hp1(2),'Color',COL2,'LineWidth',LW2,'LineStyle',LT2,...
   'Marker',MT2,'MarkerSize',MS2);
set(hp1(1),'Color',COL1,'LineWidth',LW1,'LineStyle',LT1,...
   'Marker',MT1,'MarkerSize',MS1);
if kmen1==1 & kmen2==2 & kmen3==1;
   hold on;
   hpspecial=plot(x1(:,1),x1(:,2));
   set(hpspecial,'Color',COL1,'LineWidth',LW1,'LineStyle',LT1,...
      'Marker',MT1,'MarkerSize',MS1);
   hold off;
end;
xlabel('Year','Color',COL0,'FontSize',FSlab);
ylabel('RW (mm)','Color',COL0,'FontSize',FSlab);
title('Test Plot','Color',COL0,'FontSize',FStit);
legend('Series 1','Series 2');
grid;
zoom xon;

%************  END OF DUMMY DATA SECTION 



% SOME FINAL AXIS SETTINGS

set(gca,'FontSize',FSaxis);
switch media;
case 'Paper';
   set(gca,'Color',AC1);
case 'Transp';
   set(gca,'Color',AC1);
case 'Slide';
   set(gca,'Color',AC2);
case 'Web';
   set(gca,'Color',AC2);
     
end;


% SOME FINAL FIGURE OPTIONS

switch media;
case 'Paper';
   set(gcf,'Color',FigC);
   switch orient;
   case 'Landscape';
      set(gcf,'PaperOrientation','Landscape');
      switch figlook;
      case 'Full-width';
         set(gcf,'PaperPosition',[1 1 9.5 6]);
      case 'Half-width';
         set(gcf,'PaperPosition',[1 1 9.5/2 6/2]);
      end;
   case 'Portrait';
      set(gcf,'PaperOrientation','Portrait');
      switch figlook;
      case 'Full-width';
         set(gcf,'PaperPosition',[1 1 6.5 4]);
      case 'Half-width';
         set(gcf,'PaperPosition',[1 1 3 2]);
      end;
   end;
case 'Transp';
   set(gcf,'Color',FigC);
case 'Slide';
   set(gcf,'Color',FigC,...
      'InvertHardCopy','off',...
      'PaperOrientation','Landscape',...
      'PaperPosition',[.25 .25 11.25 7.5]);
case 'Web';
   set(gcf,'Color',FigC,...
      'InvertHardcopy','off',...
   	'PaperOrientation','Landscape',...
      'PaperPosition',[.25 .25 10 7.5]);
end;
