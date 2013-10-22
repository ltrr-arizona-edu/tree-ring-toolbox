% pap2web01

kmap=questdlg('Does figure have map objects?');
kmen1=menu('Choose one',...
   'Output file for general use on web',...
   'Output file forced to slide dimension for powerpoint import');


   



% Figure background color to gray
set(gcf,'Color',[.3 .3 .3]);

% Black Line to light gray 
set(findobj(gcf,'Type','Line','Color',[0 0 0 ]),'Color',[.6 .6 .6]);

% Text to white, font 14, weight normal
set(findobj(gcf,'Type','Text'),'Color',[1 1 1],...
   'FontSize',14,...
   'FontWeight','normal');

%--- Patch edgecolors to match patch facecolors
h=findobj(gcf,'Type','patch');
npatch=size(h,1); % number of patches
if npatch>0;
   for n=1:npatch;
      hthis=h(n);
      % Any black arrow to white
      if strcmp(get(hthis,'Tag'),'Arrow');
         if all((get(hthis,'Edgecolor'))==[0 0 0]);
            set(hthis,'EdgeColor',[1 1 1]);
            set(hthis,'FaceColor',[1 1 1]);
         end;
      else; % not an arrow
         fcolor=get(h(n),'FaceColor');
         if ~strcmp(fcolor,'flat');
            set(h(n),'EdgeColor',fcolor);
         end;
      end; % end of arrow option
   end; % for n=1:npatch
end;

   
%**********  AXES PROPERTIES

ha=findobj(gcf,'Type','axes');
naxes=size(ha,1);

for n=1:naxes;
   
   h=ha(n);
   
   % Map properties
   if strcmp(kmap,'Yes');
      setm(h,'FEdgecolor',[.6 .6  .6]);  % frame color to light gray
      setm(h,'FontColor',[1 1 1]); % font color for lat/lon labels to white
   end;
   
   
   % Turn on y grid
   set(h,'YGrid','on');
   
   % Change axis background to black, and axes lines to white
   set(h,'color',[0 0 0],...
      'XColor',[1 1 1],...
      'YColor',[1 1 1],...
      'Zcolor',[1 1 1]);
   
   % Change title and axes labels to white
   set(get(h,'Title'),'Color',[1 1 1],...
      'FontSize',20,'FontWeight','normal')
   set(h,'FontSize',14);
   set(get(h,'XLabel'),'Color',[1 1 1],'FontSize',14)
   set(get(h,'YLabel'),'Color',[1 1 1],'FontSize',14)
   set(get(h,'ZLabel'),'Color',[1 1 1],'FontSize',14)
end;

if kmen1==1;
   set(gcf,'PaperOrientation','Portrait');
elseif kmen1==2; % powerpoint slide dimensions
   set(gcf,'Position',[1 29 1000*0.8  750*0.8],...
      'PaperOrientation','Portrait',...
      'PaperPosition',[1 1 10 7.5]);
end;
set(gcf,'InvertHardCopy','off');

[file1,path1]=uiputfile('*.png','Output graphics file');
pf1=[path1 file1];
eval(['print -dpng -r0 ' pf1 ';']);

% Note.  Then:
% 1 - insert into powerpoint
% 2 - save from powerpoint as png
   






