% colchg2  for re plot

set(findobj(gcf,'Color',[0 .5 0]),'Color',[0 1 0]);
set(findobj(gcf,'Type','Text'),'FontSize',12);

set(gca,'FontSize',12)

h=findobj(gcf,'Type','Text','String','T-jul_oct');
set(h,'String','T-jul-oct');
h=findobj(gcf,'Type','Text','String','P-jul_oct');
set(h,'String','P-jul-oct');

set(get(gca,'Xlabel'),'FontSize',14);
set(get(gca,'Ylabel'),'FontSize',14);
set(get(gca,'Title'),'FontSize',14);

set(findobj(gca,'Type','Line','Color',[0 1 0]),'LineWidth',2,...
   'MarkerSize',8);
set(findobj(gca,'Type','Line','Color',[0 0 1]),'LineWidth',2);

% Legend 
ha=findobj(gcf,'Type','axes');
hleg=ha(1);
set(hleg,'Color',[.2 .2 .2]);
h=findobj(hleg,'Color',[0 1 0]);
set(h,'MarkerSize',8,'LineWidth',2);
h=findobj(hleg,'Color',[0 0 1]);
set(h,'LineWidth',2);
h=findobj(hleg,'Type','Text');
set(h,'Color',[1 1 1],'FontSize',10);


set(gcf,'PaperOrientation','Landscape',...
   'PaperPosition',[.1 .1 11.25 7.5]);