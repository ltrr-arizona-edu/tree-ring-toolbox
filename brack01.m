function brack01(x,ylow,yhi,delx,dely)
% brack01: build brackets for upper and lower confidence around a point


x1 = [x-delx      x-delx    x+delx    x+delx];
ybottom = [ylow+dely   ylow      ylow      ylow+dely]; % lower
ytop    = [yhi-dely    yhi       yhi       yhi-dely]; % upper
line(x1,ybottom);
line(x1,ytop);
