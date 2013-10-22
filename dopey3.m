function respfuna(fwind,dlab,mlab,wgt,L,C,kopt,yrs)
%
% Sub-function for respfun1 that makes the page of plots
%
%************* IN
%
% fwind (1 x 1)i figure window
% dlab(2 x ?)s  label for data types (e. g., ['P';'T'])
% mlab (1 x nmos)s month labels, assumed to be in correct order; might
%   also be cell rv of season labels
% wgt (1 x 2*nmos)r  weights or correlations for the months or seas, P then T 
% L (3 x 2*nmos)L  indicates whether weights signifant at .05, .10
%     row 1 applies to .01 level.  A "1" means signif, a "0" not
%     row 2 applies to .05 level. ....
%     row 3 applies to "insignificant" correls
% C (9 x 3)r  color schemes
%   Rows 1-3 for .01, .05 and other for color verion of P
%   Rows 4-6 for .01, ..... of T
%   Rows 7-9 for ....Black and White 
% kopt -----options
%   kopt(1)  color or b/w
%      ==1 color
%      ==2 b/w
% yrs (1 x 2)i start, end year of analysis period (unadjusted N for correlations)
%
%*********** NOTES 
%
% A weights will be zero for beyond model order.  For example, if order of model
%   is 3, rows 8-9 will be zero  for the appropriate col of A
% 

% 


figure(fwind)

statusc = close(fwind);
if statusc~=1;
   close(fwind);
end

figure(fwind);

set(gcf,'RendererMode','manual');
set(gcf,'Renderer','zbuffer');


set(gcf,'DefaultLineLineWidth',2.0);
set(gcf,'DefaultTextFontWeight','bold');

if ischar(mlab);
   kmode1 = 1; % monthly model
   nmos = size(mlab,1);  % number of months in window
elseif iscell(mlab);
   kmode1 = 2; % seasonal mode
   nmos = size(mlab,2); % here mlab is a 1-row cell of season names
else;
   error('unknown mode');
end


%----------  P axes (at bottom)
ax1 = axes;
set(ax1,'Position',[.1 .09 .8 .4]);

% pull P weights
pwgt = wgt (1:nmos);
LP = L(:,1:nmos);


h1 = bar(pwgt);
x1 = get(h1,'Xdata');
y1 = get(h1,'Ydata');


%Color patches
if kopt(1)==1; % color
   c1 =C(1,:); c2 = C(2,:); c3=C(3,:);
else
   c1 =C(7,:); c2=C(8,:); c3=C(9,:);
end
CP = zeros(1,nmos,3);
n99 = sum(LP(1,:));
if n99>0;
   CP(1,LP(1,:),:) = repmat(c1,n99,1);
end
n95 = sum(LP(2,:));
if n95>0;
   CP(1,LP(2,:),:) = repmat(c2,n95,1);
end
nother = sum(LP(3,:));
if nother>0;
   CP(1,LP(3,:),:) = repmat(c3,nother,1);
end
hpatch1 = patch(x1,y1,CP);
set(hpatch1,'LineWidth',2);

% Adjust parent of patch
h1 = get(hpatch1,'parent');

set(h1,'XLim',[0 nmos+2],...
   'XTick',[1:nmos],...
   'XTickLabel',mlab,...
   'Xgrid','on',...
   'YLim',[-1 1]);

if kmode1==1;
   xlabel('Month','FontWeight','bold');
else;
   xlabel('Season','FontWeight','bold');
end

xxlim = get(gca,'XLim');
xwide = diff(xxlim);
xgo = 0.93*xwide;

text(xgo,0.70,[dlab(1,:)],'FontSize',14);
set(gca,'FontWeight','bold');

hold on;
line([0 nmos+1],[0 0],'Color',[0 0 0]);
hold off;

%Color-patch legend
xgo = xwide/100;
xdel = xwide/30;
xleg =[xgo xgo  xgo+xdel xgo+xdel];
yleg = [.75 .95 .95 .75];
hlegp1 = patch(xleg,yleg,c1);
text(xgo+xdel+xwide/100,.85,'\alpha=0.01','FontWeight','bold');
xleg =[xgo xgo xgo+xdel xgo+xdel];
yleg = [.50 .70 .70 .50];
hlegp2 = patch(xleg,yleg,c2);
text(xgo+xdel+xwide/100,.60,'\alpha=0.05');

%************** T axes (one up from bottom)
ax2 = axes;
set(ax2,'Position',[.1 .52 .8 .4]);

% pull T weights
twgt = wgt ((nmos+1):(2*nmos));
LT = L(:,(nmos+1):(2*nmos));

hb = bar(twgt);
x2 = get(hb,'Xdata');
y2 = get(hb,'Ydata');

%Color patches
if kopt(1)==1; % color
   c1 =C(4,:); c2 = C(5,:); c3=C(6,:);
else
   c1 =C(7,:); c2=C(8,:); c3=C(9,:);
end
CT = zeros(1,nmos,3);
n99 = sum(LT(1,:));
if n99>0;
   CT(1,LT(1,:),:) = repmat(c1,n99,1);
end
n95 = sum(LT(2,:));
if n95>0;
   CT(1,LT(2,:),:) = repmat(c2,n95,1);
end
nother = sum(LT(3,:));
if nother>0;
   CT(1,LT(3,:),:) = repmat(c3,nother,1);
end

% Clear bar chart and add patch
cla;
hpatch2 = patch(x2,y2,CT);
set(hpatch2,'LineWidth',2);

% Adjust parent of patch
h2 = get(hpatch2,'parent');

set(h2,'XLim',[0 nmos+2],...
   'XTick',[1:nmos],...
   'XTickLabel',' ',...
   'Xgrid','on','YLim',[-1 1]);

if kmode1==1;
   %xlabel('Month');
else;
   %xlabel('Season');
end

xxlim = get(gca,'XLim');
xwide = diff(xxlim);
xgo = 0.93*xwide;

text(xgo,0.75,[dlab(2,:)],'FontSize',14);


hold on;
line([0 nmos+1],[0 0],'color',[0 0 0]);
hold off;

set(gca,'FontWeight','bold');

%Color-patch legend
xgo = xwide/100;
xdel = xwide/30;
xleg =[xgo xgo  xgo+xdel xgo+xdel];
yleg = [.75 .95 .95 .75];
hlegt1 = patch(xleg,yleg,c1);
text(xgo+xdel+xwide/100,.85,'\alpha=0.01');
xleg =[xgo xgo xgo+xdel xgo+xdel];
yleg = [.50 .70 .70 .50];
hlegt2 = patch(xleg,yleg,c2);
text(xgo+xdel+xwide/100,.60,'\alpha=0.05');

if kmode1==1;
   txt1 = 'Monthly Correlation Analysis, ';
else;
   txt1 = 'Seasonal Correlation Analysis, ';
end
txt2 = [int2str(yrs(1)) '-' int2str(yrs(2))]; 
txt3=[txt1 txt2];
title(txt3,'FontSize',14,'FontWeight','bold');