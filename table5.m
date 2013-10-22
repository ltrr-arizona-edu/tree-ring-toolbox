function T=table5(ylog,ycal)

% compute values of flow for table 5, wr bulletin

% ylog -- 1663-1985 recon log-10 discharge
% yobs -- 1911-1985 observed log-10 discharge


% First load yhat.mat;  see ylog and yobs



ylog(:,1)=[];  % get rid of year col


L=zeros(323,9);
F=zeros (10,6);  % allocate cells of table
yr1=(1663:1985)';

L(:,1)=yr1>=1663 & yr1<=1985;
L(:,2)=yr1>=1911 & yr1<=1985;
L(:,3)=yr1>=1901 & yr1<=1975;
L(:,4)=yr1>=1891 & yr1<=1965;
L(:,5)=yr1>=1881 & yr1<=1955;

L(:,6)=yr1>=1836 & yr1<=1910;
L(:,7)=yr1>=1761 & yr1<=1835;
L(:,8)=yr1>=1686 & yr1<=1760;
L(:,9)=yr1>=1663 & yr1<=1737;


jj=[.1 .2  .5  .9];
for i=1:9
	G=ylog(L(:,i));
	F(i,1)=min(G);
	for j=1:4;
		q=jj(j);
		F(i,j+1)=quantile(G,q);
	end
	F(i,6)=max(G);
end
	
T=F;

