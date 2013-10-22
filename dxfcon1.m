function P=dxfcon1
%
% Convert .dxf file of line ends into an x,y file
% Written to get Salt and Verde streams onto a surfer map
%
% Meko 8-18-96
%
%*********************** POINT-TO FILES ***************************8
%
% file1 -- a .dxf file from autocad.  You used autocad to input
%	long/lat points, and also to input stream coordinates. You want
%	to use the coordinates to get a long-lat plot .bna plot file
%	for the stream
%
% file2 -- (out) a .mat file holding the 

% Size to hold up to 500 x,y points
a=NaN;
P=a(ones(500,1),ones(2,1));

% Get the data from file1
[file1,path1]=uigetfile('*.dxf','Input data');
pf1=[path1 file1];
fid1=fopen(pf1,'r');


%********** PART 1: READ LINES TO COUNT X,Y POINTS

na=0; % line counter
n=0; % line counter
k1=1;
while k1==1
	c=fgetl(fid1);
	na=na+1;
	if feof(fid1);
		k1=0;
	else
		nc = length(c);
		%disp(['c = ' c	]);
		%disp(['length(c) = ' int2str(length(c))]);
		if nc==4;
			if all(c=='LINE') 
				disp(['First LINE entry found at line ' int2str(na)]);
				k2=1;
				while k2==1;
					c=fgetl(fid1);
					if feof(fid1);
						k2=0;
					else
						nc=length(c);
						if nc==6;
							if all(c=='ENDSEC')
								n=n+1;
								k2=0;
							end
						else; % c~=6
								if nc==3;
									if all(c==' 10');
										n=n+1;
									end
								end
							end
						end
					end; % of if feof
				end; % of while k2
			end
		end
	end

end ; % of while k1
% Now know that have n points
ntot=n;
np=n-1; % when reach this, will do a special thing
clc
disp(['TOTAL NUMBER OF X,Y POINTS: ' int2str(ntot)]);
pause(3);



%****** PART 2: GET AND STORE X,Y POINTS

disp('STARTING PART 2')
frewind(fid1);

na=0;	% line counter
n=0; 	% line counter
k1=1;
k2=1;
while k1==1
	c=fgetl(fid1);
	na=na+1;
	nc = length(c);
	%disp(['c = ' c	]);
	%disp(['length(c) = ' int2str(length(c))]);
	if nc==4;
		if all(c=='LINE') 
			disp(['First LINE entry found at line ' int2str(na)]);
			k2=1;
			while k2==1;
				c=fgetl(fid1);
				nc=length(c);
				if nc==3;
					if all(c==' 10');
						cx=fgetl(fid1);
						x=str2num(cx);
						cdum=fgetl(fid1);
						cy=fgetl(fid1);
						y=str2num(cy);
						%[x y] 
						n=n+1;
						disp(['n = ' int2str(n)]);
						P(n,:)=[x y];
						if n==np;
							c=fgetl(fid1);
							c=fgetl(fid1);
							c=fgetl(fid1);
							cx=fgetl(fid1);
							x=str2num(cx);
							cdum=fgetl(fid1);
							cy=fgetl(fid1);
							y=str2num(cy);
							n=n+1;
							P(n,:)=[x,y];
							k2=0;
							k1=0;
						end; % of if n==np
					end; % of if all(c==' 10');
				end; % if nc==3
			end; % of while k2
		end; % of if all(c=='LINE') 
	end; % of if nc==
end; % of while k1==


fclose (fid1);

P=P(1:n,:);
