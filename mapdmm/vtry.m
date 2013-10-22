% vtry.m  
close
devicename='d:';
library = 'NOAMER';
latlim=[31 38];
lonlim=[-110 -103];
theme='bnd';
topolevel={'line'};
slnee = vmap0data(devicename,library,latlim,lonlim,theme,topolevel);
%slned = vmap0data(devicename,library,latlim,lonlim,theme,topolevel);

axesm('MapProjection','mercator','Frame','on',...
   'MapLatLimit',latlim,...
   'MapLonLimit',lonlim);
[lat,lon]=extractm(slnee,'Accurate; Definite; Coastline/Shoreline');
%displaym(slned);

plotm(lat,lon);
gridm('MLineLocation',10,'PLineLocation',10,'GLineStyle','-');
mlabel('MLabelLocation',10,'MLabelParallel','south')
plabel('PLabelLocation',10);
setm(gca,'Labelunits','dms');