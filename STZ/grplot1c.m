function grplot1c
% grplot1c: plot groups of ring-width series together on figures
% CALL: gprlot1c
% WHY: Usually run early in standardization of ring-width series.  Plots
%   help in decision of form of standarization curve to use. 
% WHO: Meko 3-20-97
%
%********************* IN ******************************************************
%
% No input arguments
%
% User prompted to point to files.  grplot1 can be run in two different modes.
%
% One mode is interactively for a single .mat ring-width file.  In this mode, user
% points to the desired .mat file, and optionally to a .mat file containing a vector
% of pointers to which series to plot and in what order.  Interactive mode allows
% zooming.
%
% Another mode is "batch", in which the user points to a file that contains filenames
% of several .rwl files;  and to another file that holds pointers for each file in
% separate rows. Batch mode does not allow zooming.
%
% In both modes, the user might optionally bypass the need for pointer files. If so,
% plots are produced for all series in the .mat ring-width storage file.
%
%
% Requirements for input files:
%  .mat ring-width file:  contains nms, yrs, X, as specified in rwlinp.m
%  .mat pointer file: contains a rv or matrix V.  Each row of V tells which ring-width
%    to plot and in what order.  V is assumed to be zero-filled to the right.  For example,
%
%    [1 5 7 8; 9 1 0 0]   means plot series 1,5,7,8 for first ring-width file, and
%       series 9 and 1 for the second ring-width file
%
%    This pointer file only needs to exist if you desire special ordering or culling
%    of ring-width series.
%   
%
%*************************** OUT *****************************************
%
% File output is optional, and consists of one or more postscript (.ps) files, one for
% each input .mat ring-width file.  The user can specify the filename in interactive 
% mode.  In batch mode, the .ps suffix is assigned to the same prefix as the input
% .mat file containing the ring-width file.  
%
% Interactive mode lets you view the figure windows as the program runs and afterwards.
% Batch mode does not save the figure windows, but overwrites them as it moves from
% ring-width file to ring-width file.
%
%****************************** UW FUNCTIONS CALLED ******************************
%
% shlmnu.m -- special menu function written by sohel anwar
%
%
%
%********************** NOTES *******************************************************8
%
% User might need to issue two commands when starting matlab session to get
% correct color and background settings for grplot1.m.  The commands are
%
% colordef none
% whitebg
%
%

%********  For colors to show correctly, and in fact for plots to show up at all
% on the axes, need to set colordef and whitebg
close all;

%*************  Use filelist or run on single mat file? *********

ButtonName=questdlg('Batch mode -- using a filelist?');
switch ButtonName
case 'No'
   filelist=0;
case 'Cancel';
   filelist=0;
case 'Yes';
   filelist=1;
   [file3,path3]=uigetfile('*.txt','Input file with list of .mat filenames');
   pf3=[path3 file3];
   if ~(exist(pf3)==2);
      error('Specified file with list of .mat filenames does not exist');
   end
end; % switch
% filelist is now set to 1 if will be using file list and
% pf3 would be the filelist file


%*************** Compute number of ring-width files to treat in this run

if filelist==0;
   nfiles=1;
else
   fid3=fopen(pf3,'r');
   nfiles=0;
   k1=1; % while control
   while k1;
      c = fgetl(fid3); 
      if ~feof(fid3);
         nfiles=nfiles+1;
      else
         k1=0;
      end
   end; % while k1
   frewind(fid3);
   
end; % filelist==0

%****************** Plot all series in each file, or use a pointer to select files **

ButtonName=questdlg('Use pointer to select ring-width series?');
switch ButtonName
case 'No'
   subset=0;
case 'Cancel';
   subset=0;
case 'Yes';
   subset=1;
   [file5,path5]=uigetfile('*.mat','Input file with selection pointer');
   pf5=[path5 file5];
   eval(['load ' pf5]); % pointer in this .mat file is assumed to be named V
   if exist('V')~=1; error(['File ' pf5 ' does not contain V']);end

   [mv,nv]=size(V);
   if (mv~=nfiles);
      error('Invalid row size for V');
   end
end; % switch
% subset is now 1 if using a selection pointer, 0 otherwise
% if subset==1, V is a rv or matrix of pointers to ring-width series

%***********************  Report current setup

disp(['Number of .mat files to process = ' int2str(nfiles)]);

% Loop over the nfile .mat files, each derived from a .rwl file

for np =1:nfiles;
   
% Get the .mat filename 
   if filelist==0; % interactively
      [flmat,path1]=uigetfile('*.mat','.MAT filename ?');
      pf1=[path1 flmat];
   else; % from file list
      c=fgetl(fid3);
      pf1=strtok(c);
      len1=length(pf1);
      fslash=findstr(pf1,'\');
      if isempty(fslash); % no path prefix in file name
         flmat=pf1;
         path1=[eval('cd') '\'];
      else
         flmat=pf1((max(fslash)+1):len1);
         path1=pf1(1:max(fslash));
      end
      % If needed, put on the .mat suffix
      if isempty(findstr(pf1,'.'));
         pf1=[pf1 '.mat'];
         flmat=[flmat '.mat'];
      end
      
   end
   
   
   % Report to screen which file working on
      disp(['   Working on .mat file: ' pf1]);
   txt1=['File # ' int2str(np) ': ' pf1];
   
   
   % Load the .mat file containing the ring-width series, names, years, etc;
   % Check that it holds the vital variables
   eval(['load ',pf1]);
   if ~all(exist('X')==1 & exist('nms')==1 & exist('yrs')==1);
      error('.mat input ring-width file does not have required variables');
   end
   
   ns=size(nms,1);  % ns is the number of ring-width series in the file
   
   if subset==1;  % Huge block applies if using vector to specify series
      
   else; % subset == 0;  Next block if not using pointer to select series
      nfigw=ceil(ns/8);     % Number of figure windows to be opened
      h=zeros(nfigw,1);
      ha=zeros(ns,1);       % Axis handles of total number of plots ns
      yht=0.1;	      % Y axis height of each plot
      k=1;
      kfp=0.01; % initial screen position (xbegin,ybegin) positions for plots
      
      % The main loop
      for i=1:nfigw,
         % Open figure windows
         h(i)=figure('Units','normal','Position',[kfp,kfp,0.85,0.85]);
         xmnx=[floor(min(yrs(k:min(ns,k+7),1))/10)*10,...
               ceil(max(yrs(k:min(ns,k+7),2))/10)*10];
         xtk=[floor(xmnx(1)/100)*100;ceil(xmnx(2)/100)*100];  % X axis limits
         xtkv=(xtk(1):100:xtk(2))';		% X tick mark vector
         ypos=0.07;	% Y position of the first plot
         for j=k:min(ns,k+7),
            ntop = min(ns,k+7); % sequential number of topmost plot
            ha(j)=axes('Position',[0.1,ypos,0.8,yht]);  % Draw the plot axes
            % Plot each individual graph
            plot((yrs(j,1):yrs(j,2))',X(yrs(j,3):yrs(j,3)+yrs(j,2)-yrs(j,1)),'b');
            %yavg=mean(X(yrs(j,3):yrs(j,3)+yrs(j,2)-yrs(j,1)));
            yrmx=max(X(yrs(j,3):yrs(j,3)+yrs(j,2)-yrs(j,1)));   % Max Y limit
            yrmn=min(X(yrs(j,3):yrs(j,3)+yrs(j,2)-yrs(j,1)));   % Min Y limit
            text('Units','normalized','Position',[0.85,0.7],'String',[int2str(j),...
                  ' - ',nms(j,:)],'Fontsize',10);    % Put the Core ID names
            ytk=[floor(yrmn);ceil(yrmx)];  % Y axis limits
            text('Units','normalized','Position',[0.1,0.85],'String',[num2str(yrmn),...
                  ' - ',num2str(yrmx)],'Fontsize',10);  % Put Y axis limit values
            % Set the X and Y axis limits and take of all axis drawings
            set(gca,'Xlim',xtk,'Ylim',ytk,'Box','off','Ytick',[]);
            set(gca,'XGrid','on'); % add verticle tick lines
            if j~=k;
               set(gca,'Xticklabel',[]);
            end;
            
            % Add title, if top plot
            if j==ntop;
               titty2 = ['  (' int2str(i) ' of ' int2str(nfigw) ')'];
               titletxt=['RINGWIDTH FILE: ' flmat];
               title([titletxt titty2]);
            end
                        
            ypos=ypos+0.1;	% starting y position for next subplot
         end
                  
         delx=1/(length(xtkv)-1);   % Y grid line interval
         xd=0;
         
         k=k+8;		% New index of plot on the next screen window
         %kfp=kfp+0.015;	% Reposition coordinates for the new figure window
         
      end; % for i =1:nfigw
      
      % ****************************** Zoom Capability
      if filelist==0;
         hzm=1;
         while hzm==1,
            % Options for Zooming
            hzm=menu('ZOOM ?','Zoom in','Zoom out','QUIT');
            if hzm==1; % zoom in
               % Build menu of possible windows to zoom
               nfigv = cellstr(num2str((1:nfigw)'));
               
               % Prompt for choosing a figure window 
               hwind=menu('Zoom on which window?',nfigv);
               figure(hwind);
               kindx=1+(hwind-1)*8;  % Row index (within nms) of first, or lowermost,
               % series in the selected zoom window
               
               % Prompt user to point to zoom region
               hhelp1 = helpdlg({'Click corners of zoom region'},...
                  'Help Message');
               pause(2);
               zpts=ginput(2);   % Get the corner points of the region to be zoomed in 
               xmx=max(zpts(:,1));
               xmn=min(zpts(:,1));
               ymx=max(zpts(:,2));
               ymn=min(zpts(:,2));
               kindx1=kindx+ceil(ymn*8);  % Minimum index of the RW series
               kindx2=kindx+floor(ymx*8); % Maximum index of the RW series
            
               xmnx=[floor(min(yrs(kindx1:min(ns,kindx2),1))/10)*10,...
                     ceil(max(yrs(kindx1:min(ns,kindx2),2))/10)*10];
               xtk=[floor((xmnx(1)+(xmnx(2)-xmnx(1))*xmn)/100)*100;...
                     ceil((xmnx(1)+(xmnx(2)-xmnx(1))*xmx)/100)*100];  % X tick vector
               xtkv=(xtk(1):100:xtk(2))';
               
               % Open a new figure window for the zoomed plot
               hzfw=figure('Units','normal','Position',[0.1,0.1,0.85,0.85]);
               axnum=(kindx2-kindx1+1);  % # of axes to be drawn
               yht=floor(80/axnum)/100;  % Y height of each plot
               ypos=0.1;		    % Initial Y position of the first plot
               for j=kindx1:kindx2,
                  haj=axes('Position',[0.1,ypos,0.8,yht]);
                  % Plot each individual graph in the zoom region
                  plot((yrs(j,1):yrs(j,2))',X(yrs(j,3):yrs(j,3)+yrs(j,2)-yrs(j,1)),'b');
                  begind=find(yrs(j,1):yrs(j,2)>=xtk(1));
                  endind=find(yrs(j,1):yrs(j,2)<=xtk(2));
                  yrmx=max(X(yrs(j,3)+begind(1)-1:yrs(j,3)+endind(length(endind))-1));
                  yrmn=min(X(yrs(j,3)+begind(1)-1:yrs(j,3)+endind(length(endind))-1));
                  text('Units','normalized','Position',[0.85,0.7],'String',...
                     [int2str(j),' - ',nms(j,:)],'Fontsize',10);	%  Core ID text
                  ytk=[floor(yrmn);ceil(yrmx)];   % Y tick mark limits
                  text('Units','normalized','Position',[0.1,0.85],'String',...
                     [num2str(yrmn),' - ',num2str(yrmx)],'Fontsize',10);
                  set(gca,'Xlim',xtk,'Ylim',ytk,'Box','off','Xtick',[],'Ytick',[]);
                  ypos=ypos+yht;
               end
               %axes('Position',[0.1,0.1,0.8,0.1]);   % Dummy axes
               % Put a Box on the whole plot
               %set(gca,'Position',[0.1,0.1,0.8,ypos-0.1],'Box','on','Ytick',[],'Xtick',[]);
               delx=1/(length(xtkv)-1);
               xd=0;
               % Draw the grid lines
               %for jj=1:length(xtkv),
                %  line('Xdata',[xd;xd],'Ydata',[0;1],'LineStyle','--','Color','k');
                %  text('Position',[xd-0.025,-0.025],'String',num2str(xtkv(jj)));
                %  xd=xd+delx;
               %end
               
            elseif hzm==2,
               close(hzfw);
               hzm=menu('ZOOM ?','Zoom in','QUIT');  % Option to zoom in other fig windows
            elseif hzm==3,
               hzm=menu('ZOOM ?','Zoom in','QUIT');  % Option to zoom in other fig windows
            end  % End of if loop
         end   % End of while loop hzm==1,
      end; % if filelist==0
   end % if subset=1
   
   
   

   
end; % for np=1:nfiles

      

disp('here');

