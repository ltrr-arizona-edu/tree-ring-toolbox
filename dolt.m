% dolt

x=rand(30,1);
y=rand(30,1);
xone  = ones(30,1)

mx=30; my=30; ncal=30; igox=1; igoy=1; 
[Ix,Iy]=indseq1(mx,my,igox,igoy,ncal,minoff);

bobs=[xone x]\y;

nsim=size(Ix,2);
B=repmat(NaN,nsim,2);

for n = 1:nsim;
   ix=Ix(:,n); iy=Iy(:,n);
   b=[xone x(ix)]\y(iy);
   B(n,:)=b';
end;

