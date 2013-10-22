function S=strX(X,tit1,nms)
[mX,nX]=size(X);

X=flipud(X);
v1 = repmat(NaN,mX,1);
v2 = repmat(NaN,mX+1,1);
X = [X v1];
X = [X ; v2'];
pcolor(X);


set(gca,'YTick',[(1:mX)+0.5],'YTickLabel',num2str(flipud([1:12]')));
set(gca,'XTick',[(1:mX)+0.5],'XTickLabel',num2str(([1:12]')));


colorbar;
title(tit1);