function uicorev(action); 
% View corev res

if nargin<1; % 
	action='initialize';
end




% Create controls
if strcmp(action,'initialize')
	close all
	start = uicontrol(gcf,...
				'Style','push',...
				'String','NULL',...
				'Position',[10 10 50 25],...
				'Visible','off',...
				'CallBack','uicorev(''startcb'')');
	pb(1) = uicontrol(gcf,...
				'Style','push',...
				'String','EPS PLOT: F2',...
				'Position',[10 50 120 25],...
				'CallBack','uicorev(''epsv'')');
	pb(2) = uicontrol(gcf,...
				'Style','push',...
				'String','QUIT',...
				'Position',[140 50 120 25],...
				'CallBack','uicorev(''closew'')');



	set(gcf,'UserData',[start pb])


	

% Define Start pushbutton callback
elseif strcmp(action,'startcb')
	ui_handles=get(gcf,'UserData');
	rb(1)=ui_handles(2);
	rb(2)=ui_handles(3);
elseif strcmp(action,'epsv')
	figure(2)
	%epsv
elseif strcmp(action,'closew')
	close all
end
