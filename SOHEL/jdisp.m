function h=jdisp(txv,pos)
%
% USAGE : h=jdisp(txv,pos)
%   This is user written function to display any text in a 
%   window for better visibility
%
%
% INPUTS
%-------
% txv	String variable to be displayed. A maximum of 2 lines
%	can be displayed. If the string is longer than 2 lines
%	it will be clipped.
% pos 	A row vector having four (4) elements specifying the
%       x,y value of the left bottom corner of the figrure
%	window and width and hight of the figure window all
% 	in normalized co-ordinates.
%
%
% OUTPUTS
%--------
% h	Display window handle
%
%
% NO USER WRITTEN FUNCTIONS NEEDED
%______________________________________________________________


if nargin==1,
  pos=[0.05 0.8 0.7 0.1];
elseif nargin==2,
  if length(pos)~=4,
    pos=[0.05 0.8 0.7 0.1];
  end
else
  pos=[0.05 0.8 0.7 0.1];
  txv=' No text input ';
end

if ~isstr(txv),
  return;
end

ls=length(txv);

figure('Numbertitle','off','Units','Normal','Position',pos,...
'Color','w','Menubar','None','Resize','off');
clf;
h=gcf;

if ls>45,
  txv1=txv(1:45);
  txv2=txv(46:ls);
  text('Units','pixels','Position',[160 35],'HorizontalAlignment','center',...
  'String',txv1,'Clipping','on','Color','k','FontWeight','Bold');
  text('Units','pixels','Position',[160 15],'HorizontalAlignment','center',...
  'String',txv2,'Clipping','on','Color','k','FontWeight','Bold');
  axis('off');
else
  text('Units','pixels','Position',[160 25],'HorizontalAlignment','center',...
  'String',txv,'Clipping','on','Color','k','FontWeight','Bold');
  axis('off');
end

% End of file 