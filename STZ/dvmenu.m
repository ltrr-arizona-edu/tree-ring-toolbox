function k = dvmenu(s0,sm);
%
% USAGE : k = dvmenu(s0,sm)
%   Generates a menu of choices for user input.
%
% K = MENU('Choose a color', A String Matrix), 
%     where String Matrix is
%	
%	    Red
%	    Blue
%	    Green
%	
% displays on the screen:
%
%	 Choose a color 
%
%	    Red
%	    Blue
%	    Green 
%
% Pressing a button on the menu window by the user returns the 
% sequence number of the choice.  On machines that support it, 
% the local menu system is used.  
%
%The maximum number of menu items is 32.
%
%
% S. Anwar, 6/17/94, revised 6/20/94 by SA.
% Copyright (c) 1994 by Control Solutions, Inc.
%________________________________________________________________

[sr,sc]=size(sm);
nriac=21;		% # of rows in a column
if sr>3*nriac,
  error('Too many items in the menu list');
  return;
end

c = computer;
display = 1;
PC = strcmp(c(1:2),'PC');
if ~strcmp(c(1:2),'PC') & ~strcmp(c(1:2),'MA')
% might be unix or VMS
	if isunix
		display = length(getenv('DISPLAY')) > 0;
	else
		display = length(getenv('DECW$DISPLAY')) > 0;
	end
end
if ~display
	disp(' ')
	disp(['----- ',s0,' -----'])
	disp(' ')
	for i=1:sr,
	    disp(['      ',int2str(i),') ',sm(i,:)])
	end
	disp(' ')
	k = input('Select a menu number: ');
	return
end

global MENU_VARIABLE;
MENU_VARIABLE = 0;
xedge = 30; %was 30;
yedge = 20; % was 35
ybord = 30 % was 30;
width = 30;
avwidth = 7; % actually 6.8886 +/- 0.4887
height = 15; % was 30
imax = 1;
maxlen = length(s0);
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
fig = figure('Position',rect,'number','off','name',' ','resize','off','Color','w',...
      'menubar','none');
set(gca,'Position',[0 0 1 1]); axis off;
% Place title ,'color','w'
t = text(mwwidth/2,mwheight-yedge/2,s0,'Horizontal','center','Vertical','top','units','pixels');
set(t,'Color','k');
set(t,'FontWeight','Bold');
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
	if PC
		waitforbuttonpress;
	else
		drawnow;
	end
end
k = MENU_VARIABLE;
delete(fig)

% End of file
