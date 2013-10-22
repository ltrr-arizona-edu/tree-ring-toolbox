function intprs2(h1,vars,lxs,ylb,dfc,N,lnt,tseg,yrY,n95,n99,indc,dy)
%
% USAGE : intprs2(h1,vars,nmX,nmY,lxs,ylb,orgn,dfc,N,lnt,...
%                 tseg,yrY,n95,n99,indc)
%   This function prints the statistical info on the figure window
%
%
% INPUTS
%-------
% h1 (1 x 1)	Figure handle
% vars 		String variable. Assumes the following values
%		'PEARSON','SPEARMAN','SIGN DEPARTURES','SIGN
%		FIRST-DIFF','HYPERGEOMETRIC'.
% lxs (1 x 1)	Length of X-segment.
% ylb 		String variable. Y-label.
% dfc		String variable matrix. Contains data type combinations.
%		like 'SKE SKE' or 'SKE RW' and so on. Also contains.
%		'Upper-sided' or 'Lower-sided' for Hypergeometric case.
% N 		String variable stating the sample size.  
% lnt (lxs x 3)	Year data for 3 top matches.
% tseg(lxs x 1)	Year segment for X-seg.
% yrY (ny x 1)	Year vector corresponding to end 
%		elements of each column of YM.
% n95 (ny x 1)	95 % confidence level.
% n99 (ny x 1)	99 % confidence level.
% indc (1 x 1)	An integer indicator.
% dy (1 x 1)	Y increment.
%
%
% NO OUTPUTS
%
%
% USER WRITTEN FUNCTIONS NEEDED
%------------------------------
% PLTEXT1.M	Places a text on the figure window
%__________________________________________________________________

ymax=0.875;	% Maximum y location in the figure window
xlf=0.15;xmd=0.45;xrt=0.75;

xlabel(['Y-year aligned with Year ',num2str(lxs),' of X-seg']);
ylabel(ylb);
pltext1(xrt,ymax,'Matches');
pltext1(xrt,ymax-dy,['1. y : ',num2str(lnt(1,1)),'-',num2str(lnt(1,1)+lxs-1)]);
pltext1(xrt,ymax-2*dy,['2. y : ',num2str(lnt(1,2)),'-',num2str(lnt(1,2)+lxs-1)]);
pltext1(xrt,ymax-3*dy,['3. y : ',num2str(lnt(1,3)),'-',num2str(lnt(1,3)+lxs-1)]);
pltext1(xlf,ymax,['X-seg :',num2str(tseg(1)),'-',num2str(tseg(lxs))]);

if strcmp(vars,'Hypergeometric '),
  if indc==1,
    dsd=dfc(2,:);
  else
    dsd=dfc(3,:);
  end
  pltext1(xmd,ymax,dsd);
end

ln9=length(n95);
lyr=length(yrY);
if ~strcmp(n95,zeros(ln9,1)),
  lnh9=round(ln9/2);
  text((yrY(lyr)+yrY(1))*0.5,n95(lnh9)*1.03,...
        '95 %','Fontsize',9);
  text((yrY(lyr)+yrY(1))*0.5,n99(lnh9)*1.03,...
        '99 %', 'Fontsize',9);
end

pltext1(xlf,ymax-dy,['Y : ',num2str(yrY(1)-lxs+1),'-',...
        num2str(yrY(lyr))]);
pltext1(xlf,ymax-2*dy, N);

% End of function
