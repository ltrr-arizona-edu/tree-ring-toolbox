function kyn=usinp(txv)
%
% USAGE : kyn=usinp(txv)
%   This is a user written routine to interactively get string input  
%   (either YES or NO) through a window. If the answer to Question 
%   'txv' is YES, kyn is set to 1 else kyn=0.
% 
% 
% INPUTS 
%-------
% txv		A text string (The question)
%
%
% OUTPUTS
%--------
% kyn (1 x 1)	kyn=1 means YES and kyn=0 means NO.
%
%
% NO USER WRITTEN FUNCTIONS NEEDED
%______________________________________________________

global mnv;
mnv=-1;

if ~isstr(txv),
  disp('Input to USINP must be a text string');
  return;
end

hf=figure('Numbertitle','off','Position',[100 300 450 100],'Color','w',...
'Menubar','None');
clf;

text('Units','pixels','Position',[160 60],'HorizontalAlignment','center',...
'String',txv,'Clipping','on','Color','k','FontWeight','Bold');
axis('off');

m=uicontrol(hf,'Style','Pushbutton','String','YES',...
  'Position',[100 10 50 30]);
set(m,'Callback',['global mnv,mnv = ',int2str(1),';']);

n=uicontrol(hf,'Style','Pushbutton','String','NO',...
  'Position',[300 10 50 30]);
set(n,'Callback',['global mnv,mnv = ',int2str(0),';']);

while mnv==-1,
  waitforbuttonpress;
end

kyn=mnv;

delete(hf)

% End of file