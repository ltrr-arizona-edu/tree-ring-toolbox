function E=gcdist(A,B)
% Great-circle distance matrix between points in A and points in B
% D Meko   9-8-93;
% rev 7-23-96 to handle decimal-degree input
%*******  INPUT ARGS
% A (mA x 2 or 4 or 6) lat and long for mA points.
%  If col dim of A is 2, cols are:
%		1 decimal longitude, negative for W
%		2 decimal latitude
%	If col dim of A is 4, cols are:
%		1 - lat deg
%		2 - lat min
%		3 - long deg
%		4 - long min
%  If col dim of A is 6, cols are:
%   	lat deg, lat min, lat sec, long deg, long min, long sec
%
% B (mB x 2 or 4 or 6) lat and long for mB  points -- see A
%***** OUTPUT ARGS
% 
% E (mA x mB) distance between points in A and points in B (km).
% 	row 1: dist from point 1 in A to points 1,2,3 ... in B
%	row 2: dist from point 2 in A "
%	etc
%******* FUNCTIONS NEEDED
%
% dms2rad -- angle from deg,min,sec to radians
% juxta1 -- offsets coords of one set of points against those 
%	of another
% dec2dms2 -- convert from decimal deg units to deg, min, sec


[mA,nA]=size(A);
[mB,nB]=size(B);


% No matter what the structure of the input A, and B, will
% want to convert to equivalent deg, min, sec

if nA==6 & nB==6; % already in dms 
elseif nA==4 & nB==4; % 4 cols, assume deg,min; expand to 6 cos
	A=[A(:,1:2) zeros(mA,1) A(:,3:4)  zeros(mA,1)];
	B=[B(:,1:2) zeros(mB,1) B(:,3:4)  zeros(mB,1)];
elseif nA==2 & nB==2; % 2-col inputs; must be decimal long and lat
	% Check that long negative; 
	L1=(all(A(:,1)<0)) & (all(B(:,1)<0)); % all longitudes negative
	if ~L1;  % a long is not negative
		error('Longs are not all negative!')
	end
	% Convert units
	A=dec2dms2(A,[1 2]); % 1 says neg west; 2 says no dialog
	B=dec2dms2(B,[1 2]);
else
	clc
	disp('A,B must have same col size, and this must be 2,4,or 6');
	error('Illegal col size for A,B')
end

D=zeros(mB,mA); % allocate
E=zeros(mA,mB);  % obsolete - see below

% Now A and B each have 6 columns
% Convert A,B to radians
Rlat1=dms2rad(A(:,1:3));
Rlong1=dms2rad(A(:,4:6));
Rlat2=dms2rad(B(:,1:3));
Rlong2=dms2rad(B(:,4:6));

U= [Rlat1 Rlong1];
V= [Rlat2 Rlong2];

%  Call function to juxtapose x,y coords (in radians) in U 
%  against those in V
Z=juxta1(U,V);


R=3968.811439/0.62137; % Denver earth radius
%R=6367.66;   % km    mean of polar and equat radii

term1=sin(Z(:,1)) .* sin(Z(:,3)) + cos(Z(:,1)) ...
.* cos(Z(:,3)) .*  cos(Z(:,2)-Z(:,4));

D1=R*acos(term1);
D(:) = D1;
E=real(D');  % Needed because matlab gives a complex number for
%    acos term if points coincide -- complex with real and imag
%    parts both of magnitude zero



