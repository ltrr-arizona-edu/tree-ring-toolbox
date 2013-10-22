function ageplot1(windnum,Y1,Y2,yrs,yrright)
% ageplot1: summary plot of tree age and time coverage by sampled rings
% CALL: ageplot1(windnum,Y1,Y2,yrs,yrright);
%
% Meko 6-22-98
%
%*********** IN
%
% windnum (1 x 1)i  window number of figure window for plot
% Y1 (? x 3)i  pith-year data
%    col 1: tree number
%    col 2: pith year, or 0 if not possible to estimate
%    col 3: code for how pith year determined:
%      1=pith present in core
%      2=pith year estimated from template for a core from this tree
%      0=not possible to estimate pith year from any core from this tree
% Y2 (? x 2)i  innermost complete ring data
%    col 1: tree number
%    col 2: year of innermost complete ring (observed) -- either by dating or ring count
% yrs (1 x 2)i  first, last year of time axis
% yrright (1 x 1)i  rightmost year at which to end bars

%****************** CHECK INPUT

if size(windnum,1)~=1 | size(windnum,2)~=1;
   error('windnum must be 1 x 1');
end

if size(Y1,1) ~= size(Y2,1);
   error('Y1 and Y2 must be same row size');
end
if size(Y1,2)~=3 | size(Y2,2)~=2;
   error('Y1 must have 3 cols, Y2 2 cols');
end
if ~all(Y1(:,1)==Y2(:,1));
   error('Col 1 of Y1 must be identical to col 1 of Y2');
end
if any(Y2(:,2)==0);
   error('Must not have missing year data in Y2');
end

[mtemp,ntemp]=size(yrs);
if mtemp~=1 | ntemp~=2;
   error('yrs must be 1 x 2');
end

[mtemp,ntemp]=size(yrright);
if mtemp~=1 | ntemp~=1;
   error('yrright must be 1 x 1');
end


ntree = size(Y1,1); % number of trees

%**************** reorder rows of Y1, Y2 so that longest lines at bottom of plot

Lzero = Y1(:,2)==0; % marks rows of Y1 that have no pith date estimated

% Check that estimated pith dates no later than first complete observed ring
yr1 = Y1(~Lzero,2);
yr2 = Y2(~Lzero,2);
if any(yr1>yr2);
   error('An estimated pith year is later than earliest complete inner ring');
end

% Build year vector with earliest of either estimated pith year or inner complete ring
yr3 = repmat(NaN,ntree,1);
yrtemp = (min([yr1'; yr2']))';
yr3(~Lzero) = yrtemp;
yr3(Lzero)=Y2(Lzero,2);

[Y3,I3]=sortrows(yr3);  % I3 is ascending ordered row index to yr3
Z1 = Y1(I3,:);   % re-ordered pith years
Z2 = Y2(I3,:);  % re-ordered innner complete observe ring

Lpith = Z1(:,2)~=0; % pointer to rows of Z1, Z3 with  pith estimate



%*********** COMPUTE PLOT PARAMETERS

yinc = 1/(ntree+2);

figure(windnum);

% Solid line for period covered by observed rings
for n = 1:ntree;
   plot([Z2(n,2) yrright],[n*yinc n*yinc],'LineWidth',1.5);
   treenum = num2str(Z2(n,1));
   text(yrright+2,n*yinc,treenum);
   if n==1;
      hold on;
   end
end

% dotted line for earlier period to estimated pith
for n = 1:ntree;
   ygo = Y3(n);
   ysp = Z2(n,2);
   plot([ygo ysp],[n*yinc n*yinc],':','LineWidth',1.5);
end

%  open circle at estimated pith, closed at observed pith
for n = 1:ntree;
   code1 = Z1(n,3);
   if code1==1;
      plot(Z1(n,2),n*yinc,'*','LineWidth',1.5);
   elseif code1==2;
      plot(Z1(n,2),n*yinc,'o','LineWidth',1.5);
   else
   end
end

xshift=diff(yrs)/30;

% Build Legend
x1=yrs(1)+5;
y1=(ntree-1)*yinc;
plot([x1 x1+32],[y1 y1],'LineWidth',1.5);
text(x1+32+xshift/2,y1,'Ring-Counted or Crossdated')

y2 = (ntree-2)*yinc;
plot([x1 x1+30],[y2 y2],':','LineWidth',1.5);
text(x1+32+xshift/2,y2,'Extrapolated');

y3 = (ntree-3)*yinc;
plot(x1+30,y3,'o','LineWidth',1.5);
text(x1+32+xshift/2,y3,'Estimated Pith Year');

y4 = (ntree-4)*yinc;
plot(x1+30,y4,'*');
text(x1+32+xshift/2,y4,'Observed Pith Year');




hold off


set(gca,'Xgrid','on','YTickLabel',[],'XLim',yrs);
