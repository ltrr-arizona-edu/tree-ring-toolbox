function P2=dxfcon2(u,v,Z1,Z2,P);
% P2=dxfcon2(u,v,Z1,Z2,P)
% 
%**************** IN ARGS ***********************************
% 
% u (1,?)r row vector of x-axis points (W to E)
% v (? x 1)r col vector of y-axis points (S to N)
% Z1 (matrix of longitudes corresponding to u,v
% Z2 (matrix of latitudes 
% P(n2 x 2)r dwg points for stream, boundary, or whatever
%
%************** OUT ARGS ***********************************
%
% P2(n3 x 2)r long-lat points for stream, boundary or whatever
%       n3==n2



[n1,m1]=size(P);
[mZ1,nZ1]=size(Z1);
[mZ2,nZ2]=size(Z2);
if mZ1~=mZ2  | nZ1~=nZ2
	error('Z1 and Z2 must be same size')
end


% Longitude axis vector -- in dwg units
	[m2,n2]=size(u);
if m2~=1,
	error('u must be row vector')
end
if n2~=nZ1;
	error('Col-size if Z1 and Z2 not equal to length of u')

end




% Latitude axis vector -- in dwg units
[m3,n3]=size(v);
if n3~=1; 
	error('v must be col vector')
end
if m3~=mZ1;
	error('Row-size of Z1 and Z2 not equal to length of v')
end


P2(:,1) =-1.0* interp2(u,v,Z1,P(:,1),P(:,2));
P2(:,2) = interp2(u,v,flipud(Z2),P(:,1),P(:,2));

% Save the ascii long-lat file
[file2,path2]=uiputfile('*.dat','Save ascii dec long-lat file as');
pf2=[path2 file2];
fid2=fopen(pf2,'w');
fprintf(fid2,'%10.4f  %10.4f\n',P2');



