function [lnt,lny,h1] = proces1(r,YM,yrY,nmX,nmY,tseg,xseg,orgn,rscl,tsm,kps,dfc,pb)
%
% USAGE : [lnt,lny,h1] = proces1(r,YM,yrY,nmX,nmY,tseg,xseg,...
%				 orgn,rscl,tsm,kps,dfc,pb)
%   A function to post-process the statistical results
%   for Pearson/Spearman, Sign, and Hypergeometric tests.
%
%
% INPUTS
%-------
% r (ny x 1)	A vector containing statistical info 
%		from Pearson/Spearman, Sign, or 
%		Hypergeometric tests
% YM (my x ny)	Time series matrix built out of vector y
%		by shifting its index by one
% yrY (ny x 1)	Year vector corresponding to end 
% nmX		String variable. Contains X-Core ID.
% nmY		String variable. Contains Y-Core ID.
% tseg(lxs x 1)	Year segment for X-seg.
% xseg(lxs x 1) X segment.
% orgn		String variable. Assumes either 'TRANSFORMED
%		SERIES' or 'ORIGINAL SERIES'.
% tsm (1 x 1)	Indicates the type of test to be performed
% kps (? x 1)	Vector containing different statistical data
% dfc		String variable matrix. Contains data type combinations.
%		like 'SKE SKE' or 'SKE RW' and so on. Also contains.
%		'Upper-sided' or 'Lower-sided' for Hypergeometric case.
% pb (1 x ny) 	Computed probability of observing m successes in 
%   		N trials
%
%
% OUTPUTS
%--------
% lnt (lxs x 3)	Year data for 3 top matches.
% lny (my x 3)	Matrix containing 3 column vectors in YM
%		corresponding to 3 best matches
% h2 (1 x 1) 	Figure handle
%
%
% USER WRITTEN FUNCTIONS NEEDED
%------------------------------
% INTPRS1.M	An intermediate post-processing function
% INTPRS2.M	An intermediate post-processing function
% INTPRS3.M	An intermediate post-processing function
% PLTEXT1.M	Places a text on the figure window
% SLHMNU.M	A modified menu function
%________________________________________________________

lxs=length(xseg);
ymax=0.875;
xlf=0.15;xmd=0.45;xrt=0.75;
nkpc=4;
 
if tsm==1,
  if kps(1)==1,
    vars='Pearson ';
  elseif kps(1)==2,
    vars='Spearman ';
  end
  nkr=length(r);
  n95=kps(3)*ones(nkr,1);
  n99=kps(5)*ones(nkr,1);
  ylb='R';
  N=['N = ',num2str(lxs)];
elseif tsm==2,
  kp=kps(1,1);
  if kp==1,
    vars='Sign-Departures ';
  elseif kp==2,
    vars='Sign-First Diff ';
  end
  nkl=length(kps);
  n95=kps(2:nkl,2);
  n99=kps(2:nkl,4);
  ylb='na';
  N=['Max ns = ',num2str(kps(1,2)),...
     ', Min ns = ',num2str(kps(1,3))];
elseif tsm==3,
  [nkpr,nkpc]=size(kps);
  vars='Hypergeometric ';
  ny=length(yrY);
  n95=kps(:,3);
  n99=kps(:,4);
  ylb='m';
  [mxm,mxi]=max(r(:,1));
  K=kps(mxi,2)*lxs;
  N=['N = ',num2str(kps(1,1))];
  q=['q = ',num2str(kps(1,2))];
  p=['p = ',num2str(pb(1))];
  if nkpc>4,
    n951=kps(:,7);
    n991=kps(:,8);
    [mxm,mxi]=max(r(:,2));
    K1=kps(mxi,6)*lxs;
    Nu=['N = ',num2str(kps(1,5))];
    qu=['q = ',num2str(kps(1,6))];
    pu=['p = ',num2str(pb(2))];
  end
end 
 
if rscl>=5,
  orgnx=' ';
else
  orgnx=orgn;
end

if nkpc<=4,
  dy=0.03;
else
  dy=0.05;
end

clrcfg;
if tsm==3 & nkpc>4,
  [rim,rsm,lnt,lny]=intprs1(r(:,1),YM,yrY,lxs);
  [rim1,rsm1,lnt1,lny1]=intprs1(r(:,2),YM,yrY,lxs);
  h1=figure('position',[50 10 550 430]);
    %        'paperposition',[0.25,1.5,7.75,8.3]);
  subplot(211);plot(yrY,r(:,1),'b',rim,rsm,'r*');
  if ~strcmp(n95,zeros(length(n95),1)),
    figure(h1); hold on;
    subplot(211);plot(yrY,n95,'g:',yrY,n99,'r-.');
  end

  title([vars,' : ',nmX,orgnx,' vs ',nmY,orgn,'  [',dfc(1,:),']']);

  intprs2(h1,vars,lxs,ylb,dfc,N,lnt,tseg,yrY,n95,n99,1,dy);
  pltext1(xlf,ymax-3*dy,['K=',num2str(K),', M=',num2str(lxs)]);
  pltext1(xmd,ymax-dy,q);
  pltext1(xmd,ymax-2*dy,p);
  pltext1(xmd,ymax-3*dy,['Max m = ',num2str(rsm(1))]);
  subplot(212);plot(yrY,r(:,2),'b',rim1,rsm1,'r*');
  if ~strcmp(n95,zeros(length(n95),1)),
    figure(h1); hold on;
    subplot(212);plot(yrY,n951,'g:',yrY,n991,'r-.');
  end
  intprs2(h1,vars,lxs,ylb,dfc,Nu,lnt1,tseg,yrY,n951,n991,2,dy);
  pltext1(xlf,ymax-3*dy,['K=',num2str(K1),', M=',num2str(lxs)]);
  pltext1(xmd,ymax-dy,qu);
  pltext1(xmd,ymax-2*dy,pu);
  pltext1(xmd,ymax-3*dy,['Max m = ',num2str(rsm1(1))]);

% Put a print pushbutton on the figure window
figure(h1);
psh=uicontrol(gcf,'Style','Pushbutton','String','PRINT','Units',...
  'Normalized','Position',[0.8 0.9 0.2 0.075],'Callback','print -v');

  bm=1;
  while bm < 4,
    bm=slhmnu('Plot Best Match ?','1st','2nd','3rd','Quit');
    if bm==4, 
      break; 
    end
    h2=intprs3(tsm,bm,lnt,lny,lxs,tseg,xseg,rsm,vars,ylb,orgn,rscl,nmX,nmY,dfc,kps,1);
    pltext1(xrt,ymax-0.06,['K=',num2str(K),', M=',num2str(lxs)]);
    h3=intprs3(tsm,bm,lnt1,lny1,lxs,tseg,xseg,rsm1,vars,ylb,orgn,rscl,nmX,nmY,dfc,kps,2);
    pltext1(xrt,ymax-0.06,['K=',num2str(K1),', M=',num2str(lxs)]);
  end
  close(h2);
  close(h3);
else
  [rim,rsm,lnt,lny]=intprs1(r,YM,yrY,lxs);
  h1=figure('position',[150 150 450 300]);
    %        'paperposition',[0.25,2,7.75,5.75]);
  plot(yrY,r,'b',rim,rsm,'r*');
  if ~strcmp(n95,zeros(length(n95),1)),
    figure(h1); hold on;
    plot(yrY,n95,'g:',yrY,n99,'r-.');
  end

  title([vars,' : ',nmX,orgnx,' vs ',nmY,orgn,'  [',dfc(1,:),']']);
  intprs2(h1,vars,lxs,ylb,dfc,N,lnt,tseg,yrY,n95,n99,1,dy);
  if tsm==3,
    pltext1(xlf,ymax-0.09,['K=',num2str(K),', M=',num2str(lxs)]);
    pltext1(xmd,ymax-0.03,q);
    pltext1(xmd,ymax-0.06,p);
    pltext1(xmd,ymax-0.09,['Max m = ',num2str(rsm(1))]);
  elseif tsm==2,
    pltext1(xmd,ymax,['Max na = ',num2str(rsm(1))]);
  end

% Put a print pushbutton on the figure window
figure(h1);
psh=uicontrol(gcf,'Style','Pushbutton','String','PRINT','Units',...
  'Normalized','Position',[0.8 0.9 0.2 0.075],'Callback','print -v');

  bm=1;
  while bm < 4,
    bm=slhmnu('Plot Best Match ?','1st','2nd','3rd','Quit');
    if bm==4, 
      break; 
    end
    h2=intprs3(tsm,bm,lnt,lny,lxs,tseg,xseg,rsm,vars,ylb,orgn,rscl,nmX,nmY,dfc,kps,1);
    if tsm==3,
      pltext1(xrt,ymax-0.06,['K=',num2str(K),', M=',num2str(lxs)]);
    end
  end
  close(h2);
end

% End of file