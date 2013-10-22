function recplot2(xdatin1,xdatin2,axin,textin1,textin2,color1)
% recplot1: tsp of reconstruction with shaded error bars and other info
% CALL: recplot1(xdatin1,xdatin2,axin,textin1,textin2,color1);
%
% Meko 8-18-00
%
% ************  IN
%
% xdatin1{}
%   {1} xr1 (nyr1 x 1)r  reconstruction, including the portion in calib period
%   {2} yrxr1 (nyr1 x 1)i year vector for xr1
%   {3} yrs (3 x 1)i
%         start year of long-term reconstruction (series xr1)
%         start year of calibration period
%         start year of smoothed period whos raw data is all in calib period
%   {4} B(nyrB x 2)upper and lower points for confidence band for xr1
%   {5} yrB (nyrB x 1) year vector for B
% xdatin2{}
%   {1} xa (nyr2 x 1)r  observed data for calibration period, possibly + later data
%   {2} yrxa (nyr2 x 1)i  year vector for xa
%   {3} xc (1 x 1)r  "event" threshold
%   {4} xm (1 x 1)r  "normal" level (e.g., mean or median)
%   {5} xytext (2 x 5)r x (row 1) and y (row 2) coords for right end of labels textin2{}
%   {6} xyarrgo (2 x 5)r  x (row 1) and y(row 2) coords for start of arrows of textin2{}
%   {7} xyarrsp (2 x 5)r  x and y coords for ends of the arrows
% axin{}
%   {1} xlimits(naxis x 2)r  limits for xaxis. naxis defines how many axes pairs
%         Starts at top plot and work down
%   {2} ylimits(naxis x 2)r limits for y axis.  If [], uses extremes of data
%   {3} pos1 (1 x 4)r  position for lowest axes (in relative coordinates)
%   {4} ysep (1 x 1)r  y-distance  separation from top of one axes pair to
%       bottom of axes above it
%   {5} Yticker (1 x ?) y tick spots
% textin1{}
%   {1} tit title for plot-- will go above top plot
%   {2} xlab xaxis label -- will go only for bottom plot
%   {3} ylab yaxis label
% textin2{} line labels (go go with arrows)
%   {1) long-term reconstructon (xr1)
%   {2} confidence band around long-term reconstruction (B)
%   {3} observed series (xa)
%   {4} threshold (xc)
%   {5} normal level (xm)
% color1(5 x 3)r 3-element RGB color specifications for elements in textin2{}
%
%*** NOTES
%
% You must add arrows graphically after producing figure
close all;

kmedia=menu('Choose media',...
   'Web',...
   'Paper');
if kmedia==1; 
   media='Web';
else;
   media='Paper';
end;


switch media;
case 'Web';
   LW1=1.5;
case 'Paper';
   LW1=1;
otherwise;
end;




%---  Unload arguments
xr1=xdatin1{1}; yrxr1=xdatin1{2}; % recon data and year vector
YRS=xdatin1{3};  % (3 x 1)i start year of recon, first year of calib period; first
%  year that xr1 data is entirely made up of data from calibration period
B=xdatin1{4}; % upper (col1) and lower (col 2) confidence bands for xr1
yrB=xdatin1{5}; % year vector for B

xa=xdatin2{1}; % observed data for calibration period, possibly + later data
yrxa=xdatin2{2}; % year vector for xa
xc=xdatin2{3}; %  "event" threshold
xm=xdatin2{4}; %  "normal" level (e.g., mean or median)
xytext=xdatin2{5}; % x,y points for labels
xyarrgo=xdatin2{6}; % x, y points for starts of arrows for labels
xyarrsp=xdatin2{7}; % x, y points of end of arrows
   

xlimits=axin{1}; % (naxis x 2)r  limits for xaxis. naxis defines how many axes pairs
%         Starts at top plot and work down
ylimits=axin{2}; % (naxis x 2)r limits for y axis.  If [], uses extremes of data
pos1=axin{3}; %  (1 x 4)r  position for lowest axes (in relative coordinates)
ysep=axin{4}; %  (1 x 1)r  y-distance  separation from top of one axes pair to
%       bottom of axes above it
yticker=axin{5};

tit=textin1{1}; %  title for plot-- will go above top plot
xlab=textin1{2};  % xaxis label -- will go only for bottom plot
ylab=textin1{3}; % yaxis label

linlab=textin2; % line labels to go with arrows; cell elements as follows:
%   {1) long-term reconstructon (xr1)
%   {2} confidence band around long-term reconstruction (B)
%   {3} observed series (xa)
%   {4} threshold (xc)
%   {5} normal level (xm)

% xyarrow; x,y points of start and end of arrows
% colorin; triplets of colors to lines or whatever in textin2 

clear xdatin1 xdatin2 textin1 textin2 axin


% Compute number of axes
naxes = size(xlimits,1);

% Compute narrowest possible y-axis limits
ydown = nanmin(nanmin([xr1; xa; B(:,1); xc; xm]));
yup = nanmax(nanmax([xr1; xa; B(:,2); xc; xm]));
ydiff=0.01 * (yup-ydown);
yup=yup+ydiff;
ydown=ydown-ydiff;

% Compute delta x between arrow and label for arrow
delx=0.005 * (xlimits(1,2)-xlimits(1,1));

% Allocate to store computed Y axis lim
YL1 = repmat(NaN,naxes,2);


%---- Build axes, bottom to top
for n = 1:naxes;
   j=naxes-n+1;
   eval(['ha' int2str(n) '= axes;']);
   pos=pos1;
   pos(2)=pos1(2)+(n-1)*[ysep + pos(4)];
   set(gca,'Position',pos);
   
   xflank=xlimits(j,:);
   
   
    % Build patch 
    colthis=color1(2,:);
   L=yrB>=xflank(1) & yrB<=xflank(2);
   if any(L);
      x1=B(L,1);
      x2=flipud(B(L,2));
      x=[x1; x2];
      yr=yrB(L);
      yr = [yr; flipud(yr)];
      hp1=patch(yr,x,colthis);
      set(hp1,'EdgeColor',colthis);
   end;
   hold on;

   
   % Plot reconstruction
   colthis=color1(1,:);
  
   L=yrxr1>=xflank(1) & yrxr1<=xflank(2);
   if any(L);
      x=xr1(L);
      yr=yrxr1(L);
      plot(yr,x,'color',colthis,'LineWidth',LW1);
   end;
   
   
   % Plot actual
   if strcmp(media,'Web');
      colthis=color1(1,:);
   else;
      colthis=[.5 .5 .5];
   end;
   L=yrxa>=xflank(1) & yrxa<=xflank(2);
   if any(L);
      x=xa(L);
      yr=yrxa(L);
      plot(yr,x,'color',colthis,'LineWidth',1.5);
   end;
  
    % Plot threshold & normal level
   colthis=color1(4,:);
   plot([xflank(1) xflank(2)],[xc xc],'color',colthis,'LineStyle','--');
   colthis=color1(5,:);
   plot([xflank(1) xflank(2)],[xm xm],'color',colthis);
   
       
   % Set xlimits
   set(gca,'Xlim',xflank);
   if ~isempty(ylimits);
      set(gca,'YLim',ylimits(j,:));
   end
   
   % X axis and observed line label on bottom plot
   if n==1;
      xlabel(xlab);
      ylabel(ylab,'HorizontalAlignment','Left');
      text(xytext(1,1),xytext(2,1),linlab{1},'HorizontalAlignment','right');
      arrow([xyarrgo(1,1) xyarrgo(2,1)],[xyarrsp(1,1) xyarrsp(2,1)],'TipAngle',10,...
         'Length',10);
      text(xytext(1,4),xytext(2,4),linlab{4},'HorizontalAlignment','center');
      arrow([xyarrgo(1,4) xyarrgo(2,4)],[xyarrsp(1,4) xyarrsp(2,4)],'TipAngle',10,...
         'Length',10);
      text(xytext(1,3),xytext(2,3),linlab{3},'HorizontalAlignment','right');
      arrow([xyarrgo(1,3) xyarrgo(2,3)],[xyarrsp(1,3) xyarrsp(2,3)],'TipAngle',10,...
         'Length',10);
   end
   
   % Title, and line label if top plot
   if n==naxes;
      switch media;
      case 'Paper';
      case 'Web';
         title(tit);
      otherwise;
      end;
      
      %text(xytext(1,3),xytext(2,3),linlab{3},'HorizontalAlignment','right');
            text(xytext(1,2),xytext(2,2),linlab{2},'HorizontalAlignment','right');
      arrow([xyarrgo(1,2) xyarrgo(2,2)],[xyarrsp(1,2) xyarrsp(2,2)],'TipAngle',10,...
         'Length',10);
           text(xytext(1,5),xytext(2,5),linlab{5},'HorizontalAlignment','right');
      arrow([xyarrgo(1,5) xyarrgo(2,5)],[xyarrsp(1,5) xyarrsp(2,5)],'TipAngle',10,...
         'Length',10);
   end

   
   % turn on x grid
   set(gca,'xgrid','on');
   
   YL1(n,:)=get(gca,'YLim');
         
   
   
end;
hold off;


%-- Re-set Y limits if desired
if isempty(ylimits);
   ylow=ydown;
   yhi = yup;
   yrng=[ylow yhi];
   for n=1:naxes;
      eval(['axes(ha' int2str(n) ');']);
      set(gca,'YLim',yrng);
      set(gca,'YTick',yticker);
   end;
   
   
end;
