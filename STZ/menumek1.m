function k = menumek1(s0,sm);
% menumek1: menu of items to click on to select one item
% CALL: k = menumek1(s0,sm)
%
% Pressing a button on the menu window by the user returns the 
% sequence number of the choice.  On machines that support it, 
% the local menu system is used.  
%
%*********************  IN ****************************
%
% s0 (1 x ?)s  the prompt to appear above the menu of choices
%    example:  'Choose a core'
% sm (sr x sc)s matrix of choices
%    Number of choices is sr
%
%************************ OUT *************************
%
% k (1 x 1)i   the sequence number of the selected item
%
%*********************** NOTES *************************
%
% The maximum number of menu items is 32.
%
%
% Meko 4-4-97; modeled on dvmenu.m,  by S. Anwar
%________________________________________________________________

global MENU_VARIABLE;
MENU_VARIABLE = 0;

% Size the string matrix of choices
[sr,sc]=size(sm);
nriac=27;		% # of rows desired in a column on the screen
ncols1=4; % maximum number of columns on the screen
if sr>ncols1*nriac,
  error('Too many items in the menu list');
  return;
end

xedge = 10; % margin from left side of menu window to left edge of 
    % menu boxes (was 30);
yedge = 12; % seems to have no effect (was 35)
ybord = 35; %ffset of top of menu from top of screen (was 30);
width = 20;  % seems to have no effect
avwidth = 7; % width of menu boxes; was 7, actually 6.8886 +/- 0.4887
height = 15; % height of menu boxes; was 30
imax = 1;
maxlen = 0;
for i = 1:sr,
    mx = length(sm(i,:));
    if mx > maxlen
       maxlen = mx;
       imax = i;
    end
end
twidth = 1.2*maxlen*avwidth;

% now figure out total dimensions needed so things can get placed in pixels
mwwidth = ceil(sr/nriac)*(twidth + 1.7*width) + 2*xedge;
mwheight = min((nriac+2),(sr+2))*yedge;
ss = get(0,'ScreenSize');
swidth = ss(3); sheight = ss(4);
%left = (swidth-mwwidth)/2;ceil(sr/nriac)
left = 20;
bottom = sheight-mwheight-ybord;
rect = [left bottom mwwidth mwheight];
fig = figure('Position',rect,'number','off','name',s0,'resize','off','Color','w',...
      'menubar','none');
set(gca,'Position',[0 0 1 1]); axis off;
% Place title ,'color','w'
%t = text(mwwidth/2,mwheight-yedge/2,s0,'Horizontal','center','Vertical','top','units','pixels');
%set(t,'Color','k');
%set(t,'FontWeight','Bold');
for ii=1:sr,
    i = (sr+1) - ii;
    nl = ceil(ii/nriac)*xedge+floor((ii-1)/nriac)*(width+twidth);
    km = min(nriac,sr);
    k = km - rem((ii-1),km);
    nb = (k-.5)*yedge;
    h1 = uicontrol('position',[nl nb width+twidth height]);
    set(h1,'callback',['global MENU_VARIABLE,MENU_VARIABLE = ',int2str(ii),';']);
    set(h1,'string',['  ', sm(ii,:)])
    set(h1,'HorizontalAlignment','left');
% left justify string inside button
end

while MENU_VARIABLE == 0
   waitforbuttonpress;
end

k = MENU_VARIABLE;
delete(fig)

% End of file
