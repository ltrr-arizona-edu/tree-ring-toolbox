%  figure4.m  --  extrapolations vs interpolations example
% residual core indices PDF7A and PDF8B
% 
load fig4dat

L7=z7(:,1)>=1794 & z7(:,1)<=1991;
L8=z8(:,1)>=1794 & z8(:,1)<=1991;


zz7=z7(L7,2);
zz8=z8(L8,2);
tt = length(zz7);
X=[ones(tt,1) zz7 zz8];

% Say calibration period is 1901 to 1991
[S,hmax]=mce1(X,[1794 1991],[1901 1991]);

te=S(:,3);
ti=~S(:,3); % interpolations

plot(zz7(te,:),zz8(te,:),'*r',zz7(ti,:),zz8(ti,:),'+r')
xlabel('Index, Tree 7','FontSize',14)
ylabel('Index, Tree 8','FontSize',14)
