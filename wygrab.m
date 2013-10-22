function wygrab
%
% Grab the water-year total flow from daily USGS data file
%
% D Meko 7-24-96
%
% Daily data obtained by USGS.  For each year, a line with 
% "WTR" in cols 2-4
% year  in cols 9-12
% data  in cols 67-74 (acre-ft for the water-year

[file1,path1]=uigetfile('*.day','Daily flow file');
fid1=fopen(file1,'r');

a=NaN;
x=a(ones(200,1),:);  % will store annual wy total in ac-ft
yr=a(ones(200,1),:); % year vector

n=0; % year counter
k1=1;

while k1==1;
	c=fgetl(fid1);
	if feof(fid1);
		k1=0;
	else; % not end of file
		nc=length(c);
		if nc<3; % no action, read next line
		else
			if   all(c(2:4)=='WTR');   % a water-year line
				n=n+1;
				syr = c(9:12);
				yr(n)=str2num(syr);
				i1=findstr(c,'AC-FT');
				igo=i1+5; % annual wy total begins after here
				c1=c(igo:nc);
				s=strtok(c1);
				x(n)=str2num(s);
				str1=sprintf('%4.0f  %7.0f',yr(n),x(n));
				disp(str1);
			else
			end
		end
	end
end

yr=yr(1:n);
x=x(1:n);
fclose(fid1);



% Check for missing years;  flag to screen if find any;  Put
% NaN where they are 

a=NaN;
yr1 = (min(yr):max(yr))';
nyrs=length(yr1);
Z=a(ones(nyrs,1),ones(2,1));
Z(:,1)=yr1;

yrgo=yr1(1);
i1 = yr-yrgo+1; 
Z(i1,2)=x;

L1= isnan(Z(:,2));
s1=sum(L1);
if s1>0;
	clc;
	disp('At least one internal year has missing WY value')
	f1=Z(L1,1);
	disp(' ');
	for n=1:s1;
		str1=sprintf('%4.0f',f1(n));
		disp(str1)
	end
	disp(' ');
	disp('Press any key to continue')
	pause
end

% Output the data in an ascii file
[file2,path2]=uiputfile('*.wyr','Output file for water-year total');

fid2=fopen([path2 file2],'w')
fprintf(fid2,'%4.0f %9.0f\n',Z');
fclose(fid2);

			





