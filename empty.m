function empty
% Expected rate of drop of water level in casing for given 
% casing diameter and pumping rate

Q=input('Pumping rate (gpm): ');
d = input('Casing diameter (in): ');

% Compute casing radius in ft
r = (d/2)/12;
A=pi*r^2;


% Convert pumping rate to units of cubic feet per minute
QQ = Q/14.7;


% Compute drop rate
drop = QQ /A ; % drop in ft/min

disp(['Drop rate (ft/min) = ',num2str(drop)])  