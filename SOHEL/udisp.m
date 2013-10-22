function udisp(txv)
%
% USAGE : udisp(txv)
%   Displays any text in a window for better visibility
%
%
% INPUTS 
%-------
% txv		A text string.
%
%
% NO OUTPUTS
%
% 
% NO USER WRITTEN FUNCTION NEEDED
%___________________________________________________________

global dv;
dv=-1;

if ~isstr(txv),
  error('Input to Udisp must be a text string');
end

ls=length(txv);

hf=figure('Numbertitle','off','Position',[50 300 450 100],'Color','w',...
'Menubar','None','Resize','off');
clf;

if ls>54,
  txv1=txv(1:50);
  txv2=txv(51:ls);
  text('Units','pixels','Position',[160 70],'HorizontalAlignment','center',...
  'String',txv1,'Clipping','on','Color','k','FontWeight','Bold');
  text('Units','pixels','Position',[160 50],'HorizontalAlignment','center',...
  'String',txv2,'Clipping','on','Color','k','FontWeight','Bold');
  axis('off');
else
  text('Units','pixels','Position',[160 60],'HorizontalAlignment','center',...
  'String',txv,'Clipping','on','Color','k','FontWeight','Bold');
  axis('off');
end

m=uicontrol(hf,'Style','Pushbutton','String','OK',...
  'Position',[200 10 50 30]);
set(m,'Callback',['global dv,dv = ',int2str(1),';']);

while dv==-1,
  waitforbuttonpress;
end

delete(hf)

% Reset default figure color to default befor leaving
set(0,'DefaultFigureColor','default');

% End of file