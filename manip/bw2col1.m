% bw2color

% Turn on y grid
set(gca,'YGrid','on');

% Change axis background to black, and axes lines to white
set(gca,'color',[0 0 0],...
   'XColor',[1 1 1],...
   'YColor',[1 1 1],...
   'Zcolor',[1 1 1]);

% Change title and axes labels to white
set(get(gca,'Title'),'Color',[1 1 1])
set(get(gca,'XLabel'),'Color',[1 1 1])
set(get(gca,'YLabel'),'Color',[1 1 1])
set(get(gca,'ZLabel'),'Color',[1 1 1])

% Black horizontal line to light gray
hline1=findobj(gca,'Type','Line','Color',[0 0 0]);
set(hline1,'Color',[.3 .3 .3]);


% Patch edgecolors to match patch colors
hp100=findobj('Type','patch','Edgecolor',[0 0 0],'FaceColor',[1 0 0]);
set(hp100,'EdgeColor',[1 0 0]);
hp177=findobj('Type','patch','Edgecolor',[0 0 0],'FaceColor',[1 .7 .7]);
set(hp100,'EdgeColor',[1 .7 .7]);
hp111=findobj('Type','patch','Edgecolor',[0 0 0],'FaceColor',[1 1 1]);
set(hp100,'EdgeColor',[1 1 1 ]);
hp001=findobj('Type','patch','Edgecolor',[0 0 0],'FaceColor',[0 0 1]);
set(hp001,'EdgeColor',[0 0 1]);
hp771=findobj('Type','patch','Edgecolor',[0 0 0],'FaceColor',[.7 .7 1]);
set(hp100,'EdgeColor',[.7 .7 1]);



% Any black text to white
htext1=findobj(gca,'Type','Text','Color',[0 0 0]);
set(htext1,'Color',[1 1 1]);

% Change figure background color to gray
set(gcf,'Color',[.3 .3 .3]);