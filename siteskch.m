function siteskch(nunits)
% siteskch:  interactive building of sketch of site collection from hand drawn sketch
% CALL: siteskch
%
% Meko 8-16-98
%
%****************** IN
%
% nunits (1 x 2)i  number of map units in x and y directions


close all;


% compute ratio of y units to x unit
rat1 = nunits(2)/nunits(1);
rat2 = 1/rat1;

% Make sure 
dely = 400;
delx = rat2 * dely;

if delx>550;
   delx = 550;
   dely = delx * rat1;
end

figpos=[50 50 delx dely];

f1=figure(1);
set(f1,'Position',figpos);


% Build axes
a1 = axes('Position',[0.1 0.1 .85 .85]);
set(gca,'XLim',[0 nunits(1)]);
set(gca,'YLim',[0 nunits(2)]);

set(gca,'XTick',0:nunits(1));
set(gca,'YTick',0:nunits(2));
set(gca,'DataAspectRatio',[1 1 1]);
grid;



