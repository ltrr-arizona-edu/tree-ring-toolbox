% Reshape a matlab array


% assume want output array X1 with   
%  nl variables per col
%  l  rows per year

nl=15;
l=17;


% assume X2 loaded


[m,n]=size(X2);

n2 = nl * l;
ndummy= n2 - n;

X3=[X2  zeros(m,ndummy)];

X4=X3';
A1=zeros(nl,l*m);
A1(:) = X4;
X1=A1';


