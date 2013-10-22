% colchg1

% Figure background color to gray
set(gcf,'Color',[.3 .3 .3]);

% Line to light gray (only lines are the horizontal mean lines
set(findobj(gcf,'Type','Line'),'Color',[.6 .6 .6]);

% Text to white, font 14, weight normal
set(findobj(gcf,'Type','Text'),'Color',[1 1 1],...
   'FontSize',14,...
   'FontWeight','normal');

%--- Patch edgecolors to match patch facecolors
h=findobj(gcf,'Type','patch');
npatch=size(h,1); % number of patches
if npatch>0;
   for n=1:npatch;
      fcolor=get(h(n),'FaceColor');
      if ~strcmp(fcolor,'flat');
          set(h(n),'EdgeColor',fcolor);
      end;
   end;
end;

   
%**********  AXES PROPERTIES

ha=findobj(gcf,'Type','axes');
naxes=size(ha,1);

for n=1:naxes;
   
   h=ha(n);
   
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

set(gcf,'Position',[1 29 1000*0.8  750*0.8],...
   'PaperOrientation','Portrait',...
   'PaperPosition',[1 1 10 7.5]);
set(gcf,'InvertHardCopy','off');

[file1,path1]=uiputfile('*.png','Output graphics file');
pf1=[path1 file1];
eval(['print -dpng -r0 ' pf1 ';']);

% Note.  Then:
% 1 - insert into powerpoint
% 2 - save from powerpoint as png
   






