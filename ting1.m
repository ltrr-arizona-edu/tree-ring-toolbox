function C=ting1(x,y)
% form a contingency table of discrete values 1,2,3,4,5

% Sample will be across top of C, master down side
mx=length(x);
C=zeros(5,5); % preallocate

A=[ones(mx,5) 2*ones(mx,5) 3*ones(mx,5)  4*ones(mx,5) 5*ones(mx,5)];
b=ones(mx,1);
b=[b 2*b 3*b 4*b 5*b];
B=[b b b b b ];

X=x(:,ones(25,1));
Y=y(:,ones(25,1));
L1=X==A;
L2=Y==B;
L3=L1 & L2;

S=sum(L3);
C(:)=S;
C(:)=flipud(C);
