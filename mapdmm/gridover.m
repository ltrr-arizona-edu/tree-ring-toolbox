function gridover(delta)
% gridover:  overlay transparency of evenly spaced grid
% CALL:  gridover;
%
%************ INPUT
%
% delta (1 x 1)r  spacing of grid lines (cm)

close all;

ha=axes;
set(ha,'DataAspectRatioMode','manual');
set(ha,'Units','centimeters');

set(ha,'Position',[1 1 15 15]);
set(ha,'XLim',[0 15],'YLim',[0 15],...
   'Xtick',[0:15],'YTick',[0:15]);

orient tall;

grid;
