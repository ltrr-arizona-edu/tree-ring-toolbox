function C=ting1(x,y)
% form a contingency table of discrete values 1,2,3

% Sample will be across top of C, master down side
mx=length(x);
C=zeros(3,3); % preallocate

A=[ones(mx,3) 2*ones(mx,3) 3*ones(mx,3)  ];
b=ones(mx,1);
b=[b 2*b 3*b ];
B=[b b b  ];

X=x(:,ones(9,1));
Y=y(:,ones(9,1));
L1=X==A;
L2=Y==B;
L3=L1 & L2;

S=sum(L3);
C(:)=S;
C(:)=flipud(C);
