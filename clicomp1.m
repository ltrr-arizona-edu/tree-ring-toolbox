function clicomp1
%
% Compare monthly climate data for a particular year at 2-4 stations
% by listing the 12 months data for the year on the screen.  Might
% need to change formats for data sets.  Written to allow easy comparison
% of Various Bisbee ppt series
%
% Meko 1-15-97

nstns=input('Number of stations to compare (2-5):  ');
a=NaN;
A=a(ones(nstns,1),ones(13,1));
b1=blanks(8);
B1=b1(ones(nstns,1),:);

year=input('Year to compare: '  );
A(:,1)=year(ones(nstns,1),:);
flag1=zeros(nstns,1); % flag that is 1 if no row for this year in a file

for n=1:nstns
	clear Z;
	% Get a file, data assumed to be in Z
	[file1,path1]=uigetfile('*.mat',['Series ' int2str(n)]);
	pf1=[path1 file1];
	eval(['load ' pf1]); % data assumed to be in Z
	if ~exist('Z')
		error(['No Z in file ' int2str(n)]);
	end

	% Store file name, padded to 8 chars if needed
	fln=strtok(file1,'.');
	nc1=length(fln);
	nmore=8-nc1;
	if nc1<8;
		badd=blanks(nmore);
		fln=[fln badd];
	else
		fln=fln(1:8);
	end
	B1(n,:)=fln;

	% Get the year of data
	yr=Z(:,1);
	L1=yr==year;
	sum1=sum(L1);
	if sum1>0;
		A(n,2:13)=Z(L1,2:13);
	else; % no row for this year in this file 
		A(n,2:13)=a(:,ones(12,1));
		flag1(n)=1;
	end
end

% Print results to screen
fmt1a='%4.0f%5.2f%5.2f%5.2f%5.2f%5.2f%5.2f';
fmt1b='%5.2f%5.2f%5.2f%5.2f%5.2f%5.2f\n';
fmt1=['%s' fmt1a fmt1b];

for n=1:nstns;
	if flag1(n)==0;
		z=A(n,:);
		b1=B1(n,:);
		str1=sprintf(fmt1,b1,z');
		disp(str1);
	else
	end
end


