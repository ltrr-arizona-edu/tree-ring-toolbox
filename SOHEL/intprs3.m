function h2=intprs3(tsm,bm,lnt,lny,lxs,tseg,xseg,rsm,vars,ylb,...
                    orgn,rscl,nmX,nmY,dfc,kps,indc)
%
% USAGE : h2=intprs3(tsm,bm,lnt,lny,lxs,tseg,xseg,rsm,vars,ylb,...
%                    orgn,nmX,nmY,dfc,kps,indc)
%   Plots and places statistical info on the figure window for best 
%   matches
% 
%
% INPUTS
%-------
% tsm (1 x 1)	Indicates the type of test to be performed
% bm (1 x 1)	Indicates which best match should be plotted
% lnt (lxs x 3)	Year data for 3 top matches.
% lny (my x 3)	Matrix containing 3 column vectors in YM
%		corresponding to 3 best matches
% lxs (1 x 1)	Length of X-segment.
% tseg(lxs x 1)	Year segment for X-seg.
% xseg(lxs x 1)	X-segment.
% rsm (3 x 1)	Vector containing 3 best values of test vector
% vars 		String variable. Assumes the following values
%		'PEARSON','SPEARMAN','SIGN DEPARTURES','SIGN
%		FIRST-DIFF','HYPERGEOMETRIC'.
% ylb 		String variable. Y-label.
% orgn		String variable. Assumes either 'TRANSFORMED
%		SERIES' or 'ORIGINAL SERIES'.
% nmX		String variable. Contains X-Core ID.
% nmY		String variable. Contains Y-Core ID.
% dfc		String variable matrix. Contains data type combinations.
%		like 'SKE SKE' or 'SKE RW' and so on. Also contains.
%		'Upper-sided' or 'Lower-sided' for Hypergeometric case.
% kps (? x 1)	Vector containing different statistical data
% indc (1 x 1)	An integer indicator
%
%
% OUTPUTS
%--------
% h2 (1 x 1) 	Figure handle
%
%
% USER WRITTEN FUNCTIONS NEEDED
%------------------------------
% PLTEXT1.M	Places a text on the figure window
%________________________________________________________

ymax=0.875;	% Maximum y location in the figure window
xlf=0.15;xmd=0.45;xrt=0.75;

if rscl>=5,
  orgnx=' ';
else
  orgnx=orgn;
end

if indc==1,
  h2=figure('position',[50 100 450 300]);
   %         'paperposition',[0.25,2,7.75,5.75]);
else
  h2=figure('position',[20 50 450 300]);
   %         'paperposition',[0.25,2,7.75,6.36]);
end

plot(lnt(:,bm),lny(:,bm),'r',lnt(:,bm),xseg,'b:');
title(['RANK ',num2str(bm),' Match             ',vars,', ',ylb,' = ',num2str(rsm(bm))]);
xlabel(['Year (',nmY,')']);
pltext1(xlf,ymax,['Dashed : ',nmX,orgnx,', ',num2str(tseg(1)),'-',num2str(tseg(lxs))]);
pltext1(xlf,ymax-0.03,['Solid    : ',nmY,orgn,', ',num2str(lnt(1,bm)),'-',num2str(lnt(1,bm)+lxs-1)]);
pltext1(xrt,ymax,dfc(1,:));
if tsm==3,
  pltext1(xrt,ymax-0.03,['N = ',num2str(kps(1,1))]);
  if indc==1,
    dsd=dfc(2,:);
  else
    dsd=dfc(3,:);
  end
  pltext1(xlf,ymax-0.06,[dsd,' q = ',num2str(kps(1,2))]);
end

% Ask if the user wants to get a particular year index from plot
kmnu=usinp('Get a year ?');
if kmnu,
  jh1=jdisp('Please click a point on the graph window using mouse');
  figure(h2);
  v=axis;
  [getyr,getln]=ginput(1);
  hold on;
  plot(getyr,getln,'r*');
  close(jh1);
  gint=find(lnt(:,bm)==round(getyr));
  gtlbx=[' x = ',num2str(tseg(gint)),' (',num2str(xseg(gint)),')'];
  gtlby=[' y = ',num2str(lnt(gint,bm)),' (',num2str(lny(gint,bm)),')'];
  vx=(getyr-v(1))/(v(2)-v(1))+0.02;
  vy=(getln-v(3))/(v(4)-v(3));
  text('Position',[vx,vy],'String',gtlbx,'Fontsize',8,...
        'Units','Normalized');
  text('Position',[vx,vy-0.03],'String',gtlby,'Fontsize',8,...
        'Units','Normalized');
end

% Put a print pushbutton on the figure window
figure(h2);
psh=uicontrol(gcf,'Style','Pushbutton','String','PRINT','Units',...
  'Normalized','Position',[0.8 0.9 0.2 0.075],'Callback','print -v');


% End of function
 