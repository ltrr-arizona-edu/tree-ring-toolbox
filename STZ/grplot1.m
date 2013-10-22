function grplot1(selvec)
%
% USAGE : gprlot1(selvec)
% This function plots multiple series offset on the same page
% 
% The mat file has a strung-out vector of data matrix in X.
% the year vector in yrs and the ID names in nms matrix. 
%
% INPUT : SELVEC - a vector containing the selected core ID indices
%
% OUTPUT : NONE
%___________________________________________________________________



% Get the .mat filename interactively
[flmat,path1]=uigetfile('*.mat','.MAT filename ?');
pf1=[path1 flmat];

% Check if the mat file exists in the workspace
if ~exist('X'),
  eval(['load ',pf1]);
end

[ns,dum1]=size(nms);  % Store the size of nms string matrix

if nargin~=0,
  if max(selvec)>ns|min(selvec)<1,
    udisp('Selected core ID index exceeds the file index');
    return;
  else
    ns=length(selvec);
  end

  nfigw=ceil(ns/8);     % Number of figure windows to be opened
  h=zeros(nfigw,1);
  ha=zeros(ns,1);       % Axis handles of total number of plots ns
  yht=0.1;	        % Y axis height of each plot
  kki=1;
  k=selvec(kki);
  kfp=0.01;

  % The main loop 

  for i=1:nfigw,
    % Open figure windows
    h(i)=figure('Color','w','Units','normal','Position',[kfp,kfp,0.85,0.85]);
  
    kmx=selvec(min(ns,kki+7));
    xmnx=[floor(min(yrs(k:kmx,1))/10)*10,ceil(max(yrs(k:kmx,2))/10)*10];
    xtk=[floor(xmnx(1)/100)*100;ceil(xmnx(2)/100)*100];  % X axis limits
    xtkv=(xtk(1):100:xtk(2))';		% X tick mark vector

    ypos=0.07;	% Y position of the first plot

    for j=kki:min(ns,kki+7),
      ha(j)=axes('Position',[0.1,ypos,0.8,yht]);  % Draw the plot axes
      % Plot each individual graph
      plot((yrs(selvec(j),1):yrs(selvec(j),2))',...
        X(yrs(selvec(j),3):yrs(selvec(j),3)+yrs(selvec(j),2)-yrs(selvec(j),1)),'b');
      %yavg=mean(X(yrs(j,3):yrs(j,3)+yrs(j,2)-yrs(j,1)));
      yrmx=max(X(yrs(selvec(j),3):yrs(selvec(j),3)+yrs(selvec(j),2)-yrs(selvec(j),1)));   % Max Y limit
      yrmn=min(X(yrs(selvec(j),3):yrs(selvec(j),3)+yrs(selvec(j),2)-yrs(selvec(j),1)));   % Min Y limit
      text('Units','normalized','Position',[0.85,0.7],'String',[int2str(selvec(j)),...
           ' - ',nms(selvec(j),:)],'Fontsize',10);    % Put the Core ID names
      ytk=[floor(yrmn);ceil(yrmx)];  % Y axis limits
      text('Units','normalized','Position',[0.1,0.85],'String',[num2str(yrmn),...
           ' - ',num2str(yrmx)],'Fontsize',10);  % Put Y axis limit values
      % Set the X and Y axis limits and take of all axis drawings
      set(gca,'Xlim',xtk,'Ylim',ytk,'Box','off','Xtick',[],'Ytick',[]);
      ypos=ypos+0.1;	% Reposition the axes for new plots
    end
    axes('Position',[0.1,ypos,0.8,yht]);   % Dummy axis
    % Put the Box on the whole plot
    set(gca,'Position',[0.1,0.07,0.8,0.9],'Box','on','Ytick',[],'Xtick',[]);
    delx=1/(length(xtkv)-1);   % Y grid line interval
    xd=0;
    % Draw the Y-grid lines along with the axis labels
    for jj=1:length(xtkv),
      line('Xdata',[xd;xd],'Ydata',[0;1],'LineStyle','--','Color','k');
      text('Position',[xd-0.025,-0.025],'String',num2str(xtkv(jj)));
      xd=xd+delx;
    end
    text('Units','normalized','Position',[0.4,0.95],'String','RINGWIDTH SERIES',...
         'Fontsize',12);    % Title of the total plots in one screen
    kki=kki+8;		% New index of plot on the next screen window
    kfp=kfp+0.015;	% Reposition coordinates for the new figure window
  end

  % Zoom Capability

  hzm=1;
  while hzm==1,
   nfigv=num2str(1);
   % Options for Zooming
   hzm=shlmnu('ZOOM ?','Zoom in','Zoom out','QUIT');
   if hzm==1,
    for i=2:nfigw,
      nfigv=str2mat(nfigv,num2str(i));  % Number of figure window vector
    end
    % Prompt for choosing a figure window 
    hfz=slvmnu('Fig window # ?',nfigv);
    figure(hfz);
    kindx=1+(hfz-1)*8;  % Starting index of the current figure window
    zdh=jdisp('Please click the corner points to zoom in');
    pause(1);
    close(zdh);
    zpts=ginput(2);   % Get the corner points of the region to be zoomed in 
    xmx=max(zpts(:,1));
    xmn=min(zpts(:,1));
    ymx=max(zpts(:,2));
    ymn=min(zpts(:,2));
    kindx1=kindx+round(ymn*8);  % Minimum index of the RW series
    kindx2=kindx+round(ymx*8); % Maximum index of the RW series

    xmnx=[floor(min(yrs(kindx1:min(ns,kindx2),1))/10)*10,...
          ceil(max(yrs(kindx1:min(ns,kindx2),2))/10)*10];
    xtk=[floor((xmnx(1)+(xmnx(2)-xmnx(1))*xmn)/100)*100;...
          ceil((xmnx(1)+(xmnx(2)-xmnx(1))*xmx)/100)*100];  % X tick vector
    xtkv=(xtk(1):100:xtk(2))';
    % Open a new figure window for the zoomed plot
    hzfw=figure('Color','w','Units','normal','Position',[0.1,0.1,0.85,0.85]);
    axnum=(kindx2-kindx1+1);  % # of axes to be drawn
    yht=floor(80/axnum)/100;  % Y height of each plot
    ypos=0.1;		    % Initial Y position of the first plot
    for j=kindx1:kindx2,
      haj=axes('Position',[0.1,ypos,0.8,yht]);
      % Plot each individual graph in the zoom region
      plot((yrs(selvec(j),1):yrs(selvec(j),2))',...
       X(yrs(selvec(j),3):yrs(selvec(j),3)+yrs(selvec(j),2)-yrs(selvec(j),1)),'b');
      begind=find(yrs(selvec(j),1):yrs(selvec(j),2)>=xtk(1));
      endind=find(yrs(selvec(j),1):yrs(selvec(j),2)<=xtk(2));
      yrmx=max(X(yrs(selvec(j),3)+begind(1):yrs(selvec(j),3)+endind(length(endind))));
      yrmn=min(X(yrs(selvec(j),3)+begind(1):yrs(selvec(j),3)+endind(length(endind))));
      text('Units','normalized','Position',[0.85,0.7],'String',[int2str(selvec(j)),...
       ' - ',nms(selvec(j),:)],'Fontsize',10);	% Put the Core ID text
      ytk=[floor(yrmn);ceil(yrmx)];   % Y tick mark limits
      text('Units','normalized','Position',[0.1,0.85],'String',[num2str(yrmn),...
           ' - ',num2str(yrmx)],'Fontsize',10);
      set(gca,'Xlim',xtk,'Ylim',ytk,'Box','off','Xtick',[],'Ytick',[]);
      ypos=ypos+yht;
    end
    axes('Position',[0.1,0.1,0.8,0.1]);   % Dummy axes
    % Put a Box on the whole plot
    set(gca,'Position',[0.1,0.1,0.8,ypos-0.1],'Box','on','Ytick',[],'Xtick',[]);
    delx=1/(length(xtkv)-1);
    xd=0;
    % Draw the grid lines
    for jj=1:length(xtkv),
      line('Xdata',[xd;xd],'Ydata',[0;1],'LineStyle','--','Color','k');
      text('Position',[xd-0.025,-0.025],'String',num2str(xtkv(jj)));
      xd=xd+delx;
    end
 
   elseif hzm==2,
    close(hzfw);
    hzm=shlmnu('ZOOM ?','Zoom in','QUIT');  % Option to zoom in other fig windows
   elseif hzm==3,
    hzm=shlmnu('ZOOM ?','Zoom in','QUIT');  % Option to zoom in other fig windows
   end  % End of if loop
  end   % End of while loop

else

  nfigw=ceil(ns/8);     % Number of figure windows to be opened
  h=zeros(nfigw,1);
  ha=zeros(ns,1);       % Axis handles of total number of plots ns
  yht=0.1;	      % Y axis height of each plot
  k=1;
  kfp=0.01;

  % The main loop 

  for i=1:nfigw,
    % Open figure windows
    h(i)=figure('Color','w','Units','normal','Position',[kfp,kfp,0.85,0.85]);
  
    xmnx=[floor(min(yrs(k:min(ns,k+7),1))/10)*10,ceil(max(yrs(k:min(ns,k+7),2))/10)*10];
    xtk=[floor(xmnx(1)/100)*100;ceil(xmnx(2)/100)*100];  % X axis limits
    xtkv=(xtk(1):100:xtk(2))';		% X tick mark vector

    ypos=0.07;	% Y position of the first plot

    for j=k:min(ns,k+7),
      ha(j)=axes('Position',[0.1,ypos,0.8,yht]);  % Draw the plot axes
      % Plot each individual graph
      plot((yrs(j,1):yrs(j,2))',X(yrs(j,3):yrs(j,3)+yrs(j,2)-yrs(j,1)),'b');
      %yavg=mean(X(yrs(j,3):yrs(j,3)+yrs(j,2)-yrs(j,1)));
      yrmx=max(X(yrs(j,3):yrs(j,3)+yrs(j,2)-yrs(j,1)));   % Max Y limit
      yrmn=min(X(yrs(j,3):yrs(j,3)+yrs(j,2)-yrs(j,1)));   % Min Y limit
      text('Units','normalized','Position',[0.85,0.7],'String',[int2str(j),' - ',nms(j,:)],...
           'Fontsize',10);    % Put the Core ID names
      ytk=[floor(yrmn);ceil(yrmx)];  % Y axis limits
      text('Units','normalized','Position',[0.1,0.85],'String',[num2str(yrmn),...
           ' - ',num2str(yrmx)],'Fontsize',10);  % Put Y axis limit values
      % Set the X and Y axis limits and take of all axis drawings
      set(gca,'Xlim',xtk,'Ylim',ytk,'Box','off','Xtick',[],'Ytick',[]);
      ypos=ypos+0.1;	% Reposition the axes for new plots
    end
    axes('Position',[0.1,ypos,0.8,yht]);   % Dummy axis
    % Put the Box on the whole plot
    set(gca,'Position',[0.1,0.07,0.8,0.9],'Box','on','Ytick',[],'Xtick',[]);
    delx=1/(length(xtkv)-1);   % Y grid line interval
    xd=0;
    % Draw the Y-grid lines along with the axis labels
    for jj=1:length(xtkv),
      line('Xdata',[xd;xd],'Ydata',[0;1],'LineStyle','--','Color','k');
      text('Position',[xd-0.025,-0.025],'String',num2str(xtkv(jj)));
      xd=xd+delx;
   end
   titletxt=['RINGWIDTH FILE: ' flmat];
    text('Units','normalized','Position',[0.4,0.95],'String',titletxt,...
         'Fontsize',12);    % Title of the total plots in one screen
    k=k+8;		% New index of plot on the next screen window
    kfp=kfp+0.015;	% Reposition coordinates for the new figure window
  end

  % Zoom Capability

  hzm=1;
  while hzm==1,
   nfigv=num2str(1);
   % Options for Zooming
   hzm=shlmnu('ZOOM ?','Zoom in','Zoom out','QUIT');
   if hzm==1,
    for i=2:nfigw,
      nfigv=str2mat(nfigv,num2str(i));  % Number of figure window vector
    end
    % Prompt for choosing a figure window 
    hfz=slvmnu('Fig window # ?',nfigv);
    figure(hfz);
    kindx=1+(hfz-1)*8;  % Starting index of the current figure window
    zdh=jdisp('Please click the corner points to zoom in');
    pause(1);
    close(zdh);
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
    hzfw=figure('Color','w','Units','normal','Position',[0.1,0.1,0.85,0.85]);
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
      text('Units','normalized','Position',[0.85,0.7],'String',[int2str(j),' - ',nms(j,:)],...
           'Fontsize',10);	% Put the Core ID text
      ytk=[floor(yrmn);ceil(yrmx)];   % Y tick mark limits
      text('Units','normalized','Position',[0.1,0.85],'String',[num2str(yrmn),...
           ' - ',num2str(yrmx)],'Fontsize',10);
      set(gca,'Xlim',xtk,'Ylim',ytk,'Box','off','Xtick',[],'Ytick',[]);
      ypos=ypos+yht;
    end
    axes('Position',[0.1,0.1,0.8,0.1]);   % Dummy axes
    % Put a Box on the whole plot
    set(gca,'Position',[0.1,0.1,0.8,ypos-0.1],'Box','on','Ytick',[],'Xtick',[]);
    delx=1/(length(xtkv)-1);
    xd=0;
    % Draw the grid lines
    for jj=1:length(xtkv),
      line('Xdata',[xd;xd],'Ydata',[0;1],'LineStyle','--','Color','k');
      text('Position',[xd-0.025,-0.025],'String',num2str(xtkv(jj)));
      xd=xd+delx;
    end
 
   elseif hzm==2,
    close(hzfw);
    hzm=shlmnu('ZOOM ?','Zoom in','QUIT');  % Option to zoom in other fig windows
   elseif hzm==3,
    hzm=shlmnu('ZOOM ?','Zoom in','QUIT');  % Option to zoom in other fig windows
   end  % End of if loop
  end   % End of while loop

end

%****************** OPTIONALLY MAKE POSTSCRIPT FILE OF THE FIGURE WINDOWS
ButtonName=questdlg('Make a postscript plot file?');
   
 switch ButtonName,
    case 'No', 
        disp('You say no');
    case 'Cancel',
         disp('You say cancel')
      case 'Yes',; % you want to make a postscript file of all figure windows
         [file2,path2]=uiputfile('*.ps','Postscript output file name');
         pf2=[path2 file2];
         for i=1:nfigw;
            figure(i)
            disp(i)
            pause(3)
            eval(['print -dps -append ' pf2]);
         end
 end % switch


% Prompt for closing all the graphs
nclg=usinp('Close Graphs ?');
if nclg,
  close all;
end


% End of file


