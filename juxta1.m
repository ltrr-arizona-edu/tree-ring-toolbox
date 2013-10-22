function Z=juxta1(A,B)
% Juxtapose x,y coordinates of 2 sets of locations

% D. Meko 9-8-93

%********  INPUT ARGS
%
% A (mA x 2) x,y coords of first set of points -- might be radians
% B (mB x 2)  likewise for second set of points


%************  OUTPUT ARGS
% 
% Z (mA*mB x 4) juxtaposed combinations of A, B
%   cols 1,2 are x,y coords of points in A
%   cols 3,4 are x,y coords of points in B


%*****  NOTES
% 
% Proceeds by duping row 1 of A mB times into a col vector.  Then
%  putting duped row 2 of A below that in same vector, etc.
%  Then for cols 3 and 4, dupe all rows of B mA times and put
%  in first nA rows of Z.  Then dupe again and put as second nA 
%  rows, etc.
%
% Use:  function called by gcdist.m in computing matrix of 
%   	great circle distances


[mA,nA]=size(A);
[mB,nB]=size(B);


% Break out the indiv cols of A;  transpose and dupe
A1=(A(:,1))';  
A2=(A(:,2))';
A1a = A1(ones(mB,1),:);
A1b = zeros(mA*mB,1);
A1b(:)=A1a;

A2a = A2(ones(mB,1),:);
A2b = zeros(mA*mB,1);
A2b(:)=A2a;

A3=[A1b A2b];


% Break out indiv cols of B;  Dupe
B1=B(:,1);
B2=B(:,2);
B1a=B1(:,ones(mA,1));
B1b=zeros(mA*mB,1);
B1b(:)=B1a;

B2a=B2(:,ones(mA,1));
B2b=zeros(mA*mB,1);
B2b(:)=B2a;

B3=[B1b B2b];

Z=[A3 B3];
