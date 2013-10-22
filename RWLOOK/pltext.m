function pltext(x,y,fontsize,txt)
% adds text txt to a plot at decimal fraction x from left side and
% decimal fraction y from bottom

v=axis;  % gets [xmin xmax ymin ymax] of current axes

x1=v(1) + x * (v(2)-v(1));   % x plotting point
y1=v(3) + y * (v(4)-v(3));   % y plotting point

text(x1,y1,txt,'FontSize',fontsize)
axis;
