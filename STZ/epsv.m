function epsv
% plot EPS vs number of trees
%load saw6
global WE WI rbtE rbtI

figure(2)
[tx,dum]=size(WI);
plot((1:tx)',WI(:,1),'y-',(1:tx)',WE(:,1),'y--');
title('EXPRESSED POPULATION SIGNAL');
xlabel('Number of Trees')
ylabel('EPS')
legend('Standard','Residual')
grid


min1 = min(min([WE(:,1) WI(:,1)]));
max1 = max(max([WE(:,1) WI(:,1)]));

set(gca,'Ylim',[min1 1])

tpos1 = (tx/10);
tpos2 = tpos1 + (tx-tpos1)/10;
ypos1 = min1+ (max1-min1)/2;
ypos2 = ypos1 - (ypos1-min1)/6;
ypos3 = ypos1 - 2*(ypos1-min1)/6;
text(tpos1,ypos1,'MEAN BETWEE-TREE CORRELATION')
rE =  round(100*rbtE(1))/100;
rI =  round(100*rbtI(1))/100;
text(tpos2,ypos3,['R: ' num2str(rE)]);
text(tpos2,ypos2,['S: ' num2str(rI)]);

figure (1)

%k1=1;
%k1 = menu('Back to GUI?','Y','N')
%	if k1==1
%		window(1)
%   else
%		window(2)
%	end
%
%close (gcf)
