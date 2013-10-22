function plotmtx1(figwind,textin,datin)
% plotmtx1:  scatterplot matrix using call to plotmatrix
% CALL: plotmtx1(figwind,textin,datin);
% 
% Meko 3-6-98
%
%**************** IN 
%
% figwind (1 x 1)i  desired figure window for plots
% textin {} text input
%  {1} (1 x 3)cell  strings with (1) bigaxis title, (2) xlabel, (3) ylabel
%  {2} (1 x ?)cell strings with x variable ids
%  {3} (1 x ?)cell strings with y variable ids
% datin {} (1 x 3)cell  
%   1- tsm of x variables (plotted along bottom)
%   2- tsm of y variables (plotted down side)
%   3- (1 x 2)i start and end year of data



% Unload data
x=datin{1};  y=datin{2};
yrgo = datin{3}(1);
yrsp = datin{3}(2);
nyr = yrsp - yrgo+1;
[mtemp1,ntemp1]=size(x);
[mtemp2,ntemp2]=size(y);
if mtemp1~=mtemp2;
   error('x and y must be same row size');
end
if nyr ~= mtemp1;
   error('row size of x and y inconsistent with start, end years in datin{3}');
end


%--------Unload text input

% title and labels for big axes 
bigtit = textin{1}(1);
xid = textin{2};
yid = textin{3};

% cell strings with x axis labels (to go along bottom)
numx = size(textin{2},2);  % number of x variables
numy = size(textin{3},2); % number of y variables

%--------  Compute correlation coefficients

R = repmat(NaN,numy,numx);
RR = corrcoef([y x]);
for n = 1:numy;
   R(n,:) = RR(n,(numy+1):(numy+numx));
end


%******************** MAKE FIGURE

figure(figwind);
[H,AX,BigAx] = plotmatrix(x,y);

[m1,n1]=size(AX);
if m1~=numy | n1~=numx;
   error('AX wrong size');
end


% label series down left (y variables)
for n=1:numy;
   a = AX(n,1);
   alaby = get(a,'ylabel');
   set(alaby,'string',yid{n});
end

% label series along bottom (x variables)
for n=1:numx;
   a = AX(numy,n);
   alabx= get(a,'xlabel');
   set(alabx,'string',xid{n});
end


%---------  Set x and y limits for the subaxes
for m = 1:numy;
   for n = 1:numx;
      a=AX(m,n);
      xlim = [min(x(:,n)) max(x(:,n))];
      ylim = [min(y(:,m)) max(y(:,m))];
      dx = abs(diff(xlim)/15);
      dy = abs(diff(ylim)/15);
      xlim = [min(x(:,n))-dx   max(x(:,n))+dx];
      ylim = [min(y(:,m))-dy   max(y(:,m))+5*dy]; % allow room to annotate r
      set(a,'XLim',xlim,'YLim',ylim);
      
      % least squares line
      axes(a); % make this axis current
      lsline;
      
      % annotated correlation coef
      str1 = sprintf('%4.2f',R(m,n));
      xstr1 = xlim(1) + dx;
      ystr1 = ylim(2) -3*dy;
      text(xstr1,ystr1,['\itr\rm = ' str1]);
    end
end

% Build title for big axes
axes(BigAx);
strb1 = sprintf('%5.0f-%5.0f (%5.0f yr)',yrgo,yrsp,nyr);
title([bigtit{1} ',  ' strb1]);
ylabel(textin{1}(3));
xlabel(textin{1}(2));

