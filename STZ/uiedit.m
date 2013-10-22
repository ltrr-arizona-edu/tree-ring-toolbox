function [x,m]=uiedit(txv,def)

% This is an m function for interactive user input 
% through editable text uicontrol.

global m k x;
%global k;
%global x;

if nargin==0,
  txv=' '; 
  def=[];
end

hf=figure('Numbertitle','off','Position',[100 200 450 100],'Color','w',...
'menubar','none','resize','off');
clf;

text('Units','pixels','Position',[160 60],'HorizontalAlignment','center',...
'String',txv,'Clipping','on','Color','k','FontWeight','Bold');
axis('off');

k=0;
x=def;

m=uicontrol(hf,'Style','Edit','Position',[200 10 40 25],...
            'String',def);
%x=get(m,'string') ,'Callback',['x = get(m,''String'');k=1;']
%pause;
set(m,'Callback',['x = get(m,''String'');']);
%set(m,'ButtonDownFcn',['x = get(m,''String'');']);

while isempty(x),
  %pause;
  set(m,'Callback',['x = get(m,''String'');']);
  %set(m,'ButtonDownFcn',['x = get(m,''String'');k=2;']);
  %x=get(m,'string');
  waitforbuttonpress;
end

h1 = uicontrol('position',[350 30 40 30]);
set(h1,'callback',['close(',int2str(hf),');']);
set(h1,'string',' OK ' );
set(h1,'HorizontalAlignment','center');

% End of file

