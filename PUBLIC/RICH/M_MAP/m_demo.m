function m_demo(num);
% M_DEMO  Demonstration program showing various maps
%


% Rich Pawlowicz (rich@ocgy.ubc.ca) 7/May/1997
% (thanks to Art Newhall for putting these examples into an m-file).
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

clf;
m_proj('ortho','lat',48','long',-123');
m_coast('patch','r');
m_grid('linest','-','xticklabels',[],'yticklabels',[]);
xlabel('Orthographic Projection','visible','on');

disp('  hit return to continue');
pause

clf;
m_proj('lambert','long',[-160 -40],'lat',[30 80]);
m_coast('patch',[1 .85 .7]);
m_elev('contourf',[500:500:6000]);
m_grid;
colormap(flipud(copper));
xlabel('Conic Projection of North America with elevations','visible','on');

disp('  hit return to continue');
pause

clf;
m_proj('stereographic','lat',90,'long',30,'radius',25);
m_elev('contour',[-3500:1000:-500],'b');
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','r');
xlabel('Polar Stereographic Projection with bathymetry','visible','on');

disp('  hit return to continue');
pause

clf;
Slongs=[-100 0;-75 25;-5 45; 25 145;45 100;145 295;100 290];
Slats= [  8 80;-80  8; 8 80;-80   8; 8  80;-80   0;  0  80];
for l=1:7,
 m_proj('sinusoidal','long',Slongs(l,:),'lat',Slats(l,:));
 m_coast('patch','g');
 m_grid('fontsize',6,'xticklabels',[],'xtick',[-180:30:360],...
        'ytick',[-80:20:80],'yticklabels',[],'linest','-','color',[.9 .9 .9]);
end;
xlabel('Interrupted Sinusoidal Projection of World Oceans','visible','on');

% The multiple maps trick is useful only with this projection. In order to
% see all the maps we must undo the axis limits set by m_grid calls:

set(gca,'xlimmode','auto','ylimmode','auto');

