function hh = slddef1
% slddef1:  slide default set. Build on colordef.m template

fig=0; % sets defaults at root level. Affects all subsequent figures

whitebg(fig,[0 0 0]);

fc=[.2 .2 .2]; % dark gray figure background
set(fig,'defaultfigurecolor',fc)
set(fig,'defaultaxescolor',[0 0 0])
set(fig,'defaultaxescolororder', ...
   1-[0 0 1;0 1 0;1 0 0;0 1 1;1 0 1;1 1 0;.25 .25 .25]) % ymcbgrw
cmap = 'defaultfigurecolormap';
set(fig,cmap,jet(64))
set(fig,'defaultsurfaceedgecolor',[0 0 0]);

set(fig,'defaultfigurepaperorientation','Portrait');
set(fig,'defaultfigurepaperposition',[0 0 11.25 7.5]);
set(fig,'defaulttextfontsize',15);

