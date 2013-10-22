function h=jdisp(txv)
%
% USAGE : h=jdisp(txv)
%   This is user written function to display any text in a 
%   window for better visibility
% I (dmm) renamed this from c:\mlb\jdisp.m 2-2-96 because sohel had a
% later 2-argument version elsewhere
%
% INPUTS
%-------
% txv	String variable to be displayed. A maximum of 2 lines
%	can be displayed. If the string is longer than 2 lines
%	it will be clipped.
%
%
% OUTPUTS
%--------
% h	Display window handle
%
%
% NO USER WRITTEN FUNCTIONS NEEDED
%______________________________________________________________

if ~isstr(txv),
  return;
end

ls=length(txv);

figure('Numbertitle','off','Position',[10 400 450 60],'Color','w',...
'Menubar','None','Resize','off');
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