% arrowdm

c=[1 0 0];
ctxt=[0 0 0];
halign='right'; % text alignment
valign='bottom';
txt='This line';
fs1=12;

harr=arrow;
xdata=get(harr,'XData');
ydata=get(harr,'YData');
hpatch=patch(xdata,ydata,c);
set(hpatch,'FaceColor',c,'EdgeColor',c);


htxt=gtext(txt);
set(htxt,'FontSize',fs1,'Color',ctxt,'HorizontalAlignment',halign,...
   'VerticalAlignment',valign);
xyztxt=get(htxt,'Position');