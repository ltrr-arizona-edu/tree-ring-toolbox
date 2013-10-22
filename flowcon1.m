function flowcon1
% flowcon1
%
% Convert a BOR-format monthly flow record (Lees Ferry, from Tom Ryan, for
% example) into a 13-col monthly file and a water-year total file.  
% 
% Meko 1-6-97
%
%******************** INFILES ***********************************
%
% *.dat -- holding the monthly data as obtained from BOR.  Data looks
%	like this for the Near Lees Ferry series (units KAF)
%
%457.900               1905-10-01
%382.150               1905-11-01
%244.550               1905-12-01
%226.530               1906-01-01
%244.810               1906-02-01
%617.510               1906-03-01
%1262.570              1906-04-01
%3960.080              1906-05-01
%5082.900              1906-06-01
%2758.190              1906-07-01
%1409.570              1906-08-01
%1430.130              1906-09-01
%
% Above, the October 1905 total was 457,900 acre-ft, etc
%
%******************** OUTFILES *************************************
%
% *.wy -- the water-year total, with year in col 1, value in col 2
% *.mon -- the monthly total, with water year in col 1, then 12 values
%
%*****************  NOTES **************************
%
% Must use the text info date at right to compute the month and water year
% Monthly file has NaN's to fill unused water year months in start and end yrs
% I have not yet tested this function on data with incomplete water years
% Output formats hard coded for pretty output with the KAF Lees Ferry data;
%	might need to change for other rivers


% Read input file and check that have start month is october (10) and end month
% is sept (9).  If not, compute number of monthly values needed to pad
% ends with NaN's to bring to full water years.  Also, count full
% (possibly padded) number of monthly values to check that divisible by
% 12 and consistent with start and end year
[file1,path1]=uigetfile('*.dat','Input monthly data file');
pf1=[path1 file1];
fid1=fopen(pf1,'r');

ncount=0;
k1=1;
cold=' ';
while k1;
	c=fgetl(fid1);
	if feof(fid1);
		clast=cold;
		k1=0;
	else
		ncount=ncount+1;
		if ncount==1;
			cfirst=c;
		end
		cold=c;
	end
end

cthis=cfirst;
fdash=find(cthis=='-');
yrgo = str2num(cthis((fdash(1)-4):(fdash(1)-1)));
mongo = str2num(cthis((fdash(2)-2):(fdash(2)-1)));

cthis=clast;
fdash=find(cthis=='-');
yrsp = str2num(cthis((fdash(1)-4):(fdash(1)-1)));
monsp = str2num(cthis((fdash(2)-2):(fdash(2)-1)));

ngo=0; % pad start with this many monthly NaN's
nsp=0; % pad end with etc
if mongo~=10;
	disp('Start month not October -- will pad');
	if mongo>10;
		ngo=mongo-10;
		wygo=yrgo+1;
	else; % mongo is less than 10
		nsp=2+mongo;
		wygo=yrgo;
	end
else
	wygo=yrgo+1;  % starting water year for output
end
	
if monsp~=9;
	disp('End month not Sept -- will pad');
	if monsp>9;
		ngo=9+(12-monsp);
		wysp=yrsp+1;
	else; % monsp is less than 9
		nsp=9-monsp;
		wysp=yrsp;
	end
else
	wysp=yrsp;
end

montot=ncount+ngo+nsp;
if rem(montot,12);
	error('Total number of monthly values not divisible by 12');
end

% Size the monthly matrix 
yr=(wygo:wysp)';
nyrs = wysp-wygo+1;
a=NaN;
X=a(ones(nyrs,1),ones(13,1));
X(:,1)=(wygo:wysp)';

% Size the temporary monthly vector
y=a(ones(ncount,1),:);


% Rewind and re-read the input file, monthly values before padding 
frewind(fid1);
for n = 1:ncount;
	c=fgetl(fid1);
	y(n)=str2num(strtok(c));
end

% Pad front if needed
if ngo>0;
	w=a(ones(ngo,1),:);
	y=[w; y];
end

% Pad end if needed
if nsp>0;
	w=a(ones(nsp,1),:);
	y=[y; w;];
end


% Reshape the monthly data into a 12-col matrix, organized by water year
Z=reshape(y,12,nyrs);
Z=Z';

% Compute the water year totals
z=(sum(Z'))';

% Make the output monthly matrix
X=[yr Z];

% Make the output wy totals matrix
x=[yr z];


% Make the output files -- fid2 for the monthly, fid3 for the wy total
[file2,path2]=uiputfile('*.mon','File with monthly totals, output');
pf2=[path2 file2];
fid2=fopen(pf2,'w');
[file3,path3]=uiputfile('*.wy','File with wy totals output');
pf3=[path3 file3];
fid3=fopen(pf3,'w');

fmt1='%4.0f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f\n';
fmt2='%4.0f%10.3f\n'

fprintf(fid2,fmt1,X');
fprintf(fid3,fmt2,x');

fclose all
