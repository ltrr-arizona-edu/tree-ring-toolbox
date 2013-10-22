function zoomspec (f,y1,y2,yl,yh)

% Zoom on a pdgm1.m -produced spectrum

pd1=input('Longest period in window:')
pd2=input('Shortest period in window:')

g1 = 1 ./ pd1;   % lowest freq
g2 = 1 ./ pd2;   % highest freq

L=f <= g2 & f >= g1;

plot(f(L),y1(L),'-g',f(L),yl(L),'r:',f(L),y2(L),'--w')
ylabel('Spectral Estimate')
title('Methuselah Walk, 6000 BC - AD 1979')

  
pause

np =input('How many peaks to mark?:')

for j=1:np;
	[x,y]=ginput(2);
	L1= f>= x(1) & f<=x(2);
	f1=f(L1);
	s1=y1(L1);
	yd=y1(L1) - y2(L1);
	[ydmax,imax]=max(yd);
	fmax=f1(imax);
	pmax=1/fmax;
	gtext(num2str(pmax));
end

