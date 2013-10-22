function rwwind01(figspecs,strin,datin,figwind);
% rwwind01(figspecs,datin) is a dedicated subfunction of rwlook that builds a
% print window for user to take results to microscope

id1=strin{1};
id2=strin{2};
title5=strin{3};

huse=figspecs{1};  % part of height of window to use (.95 is reasonable)
w = figspecs{2};  % width of figures, as fraction of window
d = figspecs{3}; % vertical distance between plots (gap)
p1 = figspecs{4}; % x, y coordinate of lower left of lower plot
f = figspecs{5}; % fraction height of lower window is of height of upper (e.g., 0.25)


t1=datin{1};
t2=t1;
t2(1)=[];
wx=datin{2};
wy=datin{3};
dx = datin{4};
dy = datin{5};

yreh=datin{6};
v=datin{7};
nH=datin{8};
tau95=datin{9};
tau99=datin{10};

xlist=datin{11};
ylist=datin{12};


% Compute heights of the 6 subfigures
havail = huse - 2*d; % Compute height available to axes
h5=havail/(4*f+1);  % height of figure 5
h1 = f*h5;  % height of figs 1,2,3,4
h2=h1; h3=h1; h4=h1;
h6= h5+4*h1; % for listing of ringwidths

% Compute the y position for lower edge of axes
y1=p1(2);  
y2=y1+h1;
y3=y2+h1+d/2;
y4=y3+h1;
y5=y4+h1+d;
y6=y1;

% compute x postion
x1=p1(1); x2=x1; x3=x1; x4=x1; x5=x1; x6=w+p1(1);

% Compute position vectors for axes
pos1 = [ x1 y1 w h1];
pos2 = [x2 y2 w h2];
pos3=[x3 y3 w h3];
pos4=[x4 y4 w h4];
pos5=[x5 y5 w h5];
pos6=[x6 y6 (1-w) h6];

%************  BUILD FIGURE

figure(figwind);

% Rw of master series y
axes1=axes('Position',pos1);
stem(t1,wy);
set(axes1,'XGrid','on','YGrid','off',...
   'TickLength',[0 .025]);
ylabel(id2);
set(gca,'TickDir','out');

% RW of test series x
axes2=axes('Position',pos2);
stem(t1,wx);
set(axes2,'XGrid','on','YGrid','off',...
   'TickLength',[0 .025],'XTickLabel','');
% Annotate 'RW' in upper left of plot
xlim2=get(axes2,'XLim');
ylim2=get(axes2,'YLim');
delx = diff(xlim2);
dely = diff(ylim2);
text(xlim2(1)+delx/100,ylim2(2)-dely/8,'RW');
ylabel(id1);

% scaled rw change, master series
axes3=axes('Position',pos3);
Lup = dy>0;
Ldn = dy<=0;
hp3a = stem(t2(Lup),dy(Lup));
set(hp3a,'color','b');
hold on;
hp3b = stem(t2(Ldn),dy(Ldn));
set(hp3b,'color','r');
hold off;
line(t2,zeros(length(t2),1));
set(axes3,'XGrid','on','YGrid','off',...
   'TickLength',[0 .025],'XTickLabel','');
ylabel(id2);
set(gca,'XLim',xlim2);



% scaled rw change, test series
axes4=axes('Position',pos4);
Lup = dx>0;
Ldn = dx<=0;
hp4a = stem(t2(Lup),dx(Lup));
set(hp4a,'color','b');
hold on
hp4b = stem(t2(Ldn),dx(Ldn));
set(hp4b,'color','r');
hold off;
line(t2,zeros(length(t2),1));
set(gca,'XLim',xlim2);
ylabel(id1);

set(axes4,'XGrid','on','YGrid','off',...
   'TickLength',[0 .025],'XTickLabel','');
% Annotate 'delta RW' in upper left of plot
xlim4=get(axes4,'XLim');
ylim4=get(axes4,'YLim');
delx = diff(xlim4);
dely = diff(ylim4);
text(xlim4(1)+delx/100,ylim4(2)-dely/8,'\DeltaRW/g');

axes5=axes('Position',pos5);
plot(yreh,v,yreh,tau95(ones(nH,1),:),'m--',...
   yreh,tau99(ones(nH,1),:),'m--');
title(title5);

% Window listing ring widths
axes6=axes('Position',pos6);
fmt1 = '%4.0f   %3.0f  %3.0f\n';
nvals = length(t1);
yinc = 1/(nvals+3);
header1 = [strtok(id1) ' ' strtok(id2)];
text(0.1,1+yinc,header1,'FontSize',9);
for nn = 1:nvals;
   strlist = sprintf(fmt1,t1(nn),100*xlist(nn),100*ylist(nn));
   text(0.1,(1-yinc*nn),strlist,'FontSize',9);
end

set(gca,'Visible','off');
    