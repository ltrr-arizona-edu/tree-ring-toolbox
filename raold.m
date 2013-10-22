function respfuna(fwind,dlab,mlab,wgt,L,A,C,kopt)
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
% A(8 x 2)r  autoregressive model results.  Col 1 for long period, col 2 for short
%   Rows:
%   1 start year
%   2 end year
%   3 ar model order
%   4 pct variance accounted for by autoregression
%   5-9  acf, lags 1-5 
%   10-14 2 SE of acf at lags 1-5
% C (9 x 3)r  color schemes
%   Rows 1-3 for .01, .05 and other for color verion of P
%   Rows 4-6 for .01, ..... of T
%   Rows 7-9 for ....Black and White 
% kopt -----options
%   kopt(1)  color or b/w
%      ==1 color
%      ==2 b/w
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
set(ax1,'Position',[.1 .09 .8 .3]);

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
xgo = 0.70*xwide;

text(xgo,0.70,['Correlation: ',dlab(1,:)],'FontSize',14);
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
set(ax2,'Position',[.1 .42 .8 .3]);

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
xgo = 0.70*xwide;

text(xgo,0.75,['Correlation: ',dlab(2,:)],'FontSize',14);


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



%  AR model weights, long-period model
ax3 = axes;
set(ax3,'Position',[.13 .75 .3 .2]);

hstem1 = stem([1:5]',A(5:9,1));
hold on;
hs1a = plot([1  5],[0 0],[1:5],A(10:14,1),[1:5],-A(10:14,1));
hold off;
set(hs1a,'color',[0 0 0]);
set(gca,'YLim',[-1 1],'XLim',[0 5.5],...
   'XTick',[1 2 3 4 5],...
   'XTickLabel',[]);
line([0 5.5],[0 0],'Color',[0 0 0]);   
str1 = sprintf('%4.0f-%4.0f',A(1:2,1));
title(['Tree Rings, ' str1],'FontWeight','bold');
str2a = sprintf('%3.1f',100*A(4,1));
str2 = ['AR(' int2str(A(3,1)) '): ' str2a '%'];
text(1.5,-0.8,['Model: ',str2],'FontSize',8);
text(3,0.8,'ACF and 2 SE','FontSize',8);
set(gca,'FontWeight','bold');
   

%  AR model weights, instrumental period
ax4 = axes;
set(ax4,'Position',[.53 .75 .3 .2]);
hstem2 = stem([1:5]',A(5:9,2));
hold on;
hs2a = plot([1  5],[0 0],[1:5],A(10:14,2),[1:5],-A(10:14,2));
hold off;
set(hs2a,'color',[0 0 0]);

set(gca,'YLim',[-1 1],'XLim',[0 5.5],...
   'XTick',[1 2 3 4 5],...
   'XTickLabel',[]);
line([0 5.5],[0 0],'Color',[0 0 0]);   
str1 = sprintf('%4.0f-%4.0f',A(1:2,2));
title(['Tree Rings, ' str1],'FontWeight','bold');
str2a = sprintf('%3.1f',100*A(4,2));
str2 = ['AR(' int2str(A(3,2)) '): ' str2a '%'];
text(1.5,-0.8,['Model: ',str2],'FontSize',8);
text(3,0.8,'ACF and 2 SE','FontSize',8);
set(gca,'FontWeight','bold');
str2;
