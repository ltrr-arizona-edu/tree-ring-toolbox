yrx=(1916:1990)';
len=50;
offset=5;

[yrh,H]=pullseg1(h,yrx,len,offset);
[yrb,B]=pullseg1(bh,yrx,len,offset);

[mH,nH]=size(H);

z=zeros(nH,1);

for i=1:nH;
	x=H(:,i);
	y=B(:,i);
	[T,tau]=kendtau(x,y);
	z(i)=tau;
end

plot(yrb,z);
pltext(.1,.9,[num2str(len),'-year segments']);
pltext(.1,.8,[' Offset by ',num2str(offset),' years']);
ylabel('Kendall tau');
xlabel('Ending year');