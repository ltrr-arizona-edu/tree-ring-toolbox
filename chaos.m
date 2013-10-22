hold off
n=35;
rand('normal');
x=zeros(n,20);
x(1,1:20)=rand(1,20);


for t=2:n;
	x(t,:)=a*x(t-1,:) .* (1-x(t-1,:));
end

hold on
for j=1:20
	plot(x(:,j))
end
