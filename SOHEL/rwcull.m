function [xseg,tseg] = rwcull(x,y,yx,yy,nmx,nmy,orgn,rscl)
%
% USAGE : [xseg,tseg] = rwcull(x,y,yx,yy,nmx,nmy,orgn)
%   Returns a segment of x data series xseg defined by the user
%   The user is asked to enter the segment specification either
%   through graphical input or through specification of beginning
%   and ending years of the segment.
%
%
% INPUTS
%-------
% x (nx x 1)    : x data series.
% y (ny x 1)    : y data series.
% yx (nx x 1)   : Year vector for x.
% yy (ny x 1)   : Year vector for y.
% nmx	  	: String variable. Core ID for x series.
% nmy	  	: String variable. Core ID for y series.
%
%
% OUTPUTS
%--------
% xseg(lxs x 1) : Segment of X-series chosen.
% tseg(lxs x 1) : Year vector corresponding to xseg.
%
% 
% USER WRITTEN FUNCTIONS NEEDED 
%------------------------------
% CLSGRF.M	A script .m file called by uicontrol which closes
%		all figure windows.
% CLRCFG.M	A script file to set white background
% JDISP.M	A non-interactive window display function
% PLTEXT.M	Places text inside figures.
% SHLMNU.M	A modified menu function
%_________________________________________________________________

% Set the root default colors
clrcfg;

if rscl>=5,
  orgnx=' ';
else
  orgnx=orgn;
end

% Always plot initial X and Y series in figure 1
figure(1);
clf;
subplot(211),plot(yx,x,'k'),title(['Time Series Xall : ',nmx,orgnx]);
xlabel('Tentative years'),ylabel('Xall');
subplot(212),plot(yy,y,'m'),title(['Time Series Y : ',nmy,orgn]);
xlabel('Years'),ylabel('Y');

% Plot only the X series in figure 2 (h)
h=figure('position',[150 150 450 300],...
         'paperposition',[0.25,2,7.75,5.75]);
clf;
subplot(111), plot(yx,x,'b'),title(['Time Series Xall : ',nmx,orgnx]);
xlabel('Tentative years'),ylabel('Xall');

% Keep the menu color as default
k=shlmnu('Choose X-segment','Specify years','Graph input');

if k==1,
  jhl=jdisp('Please Return to Command Window to Enter Years');
  tbeg=input('Enter the Start year : ');
  tend=input('Enter the end year : ');
  close(jhl);
   % Check for the correct start and end year
   while tbeg>tend | tbeg>max(yx), 
    tbeg=input('Please enter the correct Start year : ');
    tend=input('Please enter the  correct end year : ');
   end
else
  [t,v]=ginput(2);
  tbeg=ceil(t(1));
  tend=floor(t(2));
   % Check for the correct start and end year
   while tbeg>tend | tbeg>max(yx),   
    dh=jdisp(['Start or End year out of range.',...
    'Please Specify the correct segment.']);
    pause(3);
    figure(h);
    [t,v]=ginput(2);
    tbeg=ceil(t(1));
    tend=floor(t(2));
    close(dh); % close the jdisp window
   end
end

% If the input start year is less than minimum year in the data,
% set the start year as minimum value in data by default
if tbeg<min(yx), 
  tbeg=min(yx); 
  disp('Start year out of range : Plotting with default value');
end
% If the input end year is less than maximum year in the data,
% set the end year as maximum value in data by default
if tend>max(yx), 
  tend=max(yx); 
  disp('End year out of range : Plotting with default value');
end

tseg=(tbeg:tend)';
ibeg=find(yx==tbeg);
iend=find(yx==tend);
xseg=x(ibeg:iend);

% Always plot X segment and Y series in figure 1
figure(1);
clf;
tbegs=num2str(tbeg);tends=num2str(tend);
subplot(211),plot(tseg,xseg,'r');
title(['Time Series : ',nmx,orgnx,' Segment : ',tbegs,'-',tends]);
xlabel('Segment Year'),ylabel('Xseg');
subplot(212),plot(yy,y,'m');
title(['Time Series Y : ',nmy,orgn,' : ',num2str(yy(1)),'-',...
       num2str(yy(length(yy)))]);
xlabel('Year'),ylabel('Y');

% Put a print pushbutton on the figure window
figure(1);
psh=uicontrol(gcf,'Style','Pushbutton','String','PRINT',...
  'Position',[400 350 120 30],'Callback','print -v');

close(h);  % Close the X series figure window before leaving

% End of file

