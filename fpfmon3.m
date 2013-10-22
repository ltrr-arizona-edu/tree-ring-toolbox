function fpfmon3(Y,fid1,stn,ste,id,dtype,dunits,z,digs)
% Print monthly climate data in "SRP-appendix" format
%
% D Meko 7-6-96
%
%********************  IN ARGS  ***************************
%
% Y (mY x 13)r   year in col 1, 12 monthly data 
%		values in cols 2-13
% fid1 (?)  handle for file to write data to
% stn (1 x 70)s station name
% ste (1 x 70)s state
% id (1 x 70)s station id
% dtype (1 x 1)i  data type:
%		1=pcp
%		2=tmp
%		3=streamflow
%		4=evap
% dunits (1 x 70)s  data units
% z (1 x 3)r  long, lat, elev -- dec degrees, meters

% tit (1 x 70)s 4 title 
% digs (1 x 1)i scaling digit, typically 2 for ppt data
%    Data is multiplied by 1E"digs", then rounded to nearest integer;
%    Digs=2 gives 1E2=100 as scaling, which produces hundredths of 
%	  inches if Y is input as inches
%
%***************************  OUT ARGS -- NONE
%
%********************  NOTES  **********************
% 
% Scales and rounds so that printed data has no decimal points
% For example, pcp might be printed in hundredths of inches, so that
%		1.29 inches prints as 129
% Year and month headings hard-coded
% Numeric field is 5, allowing for max of 99.99 inches (9999) with
%  one blank spacing between months
% Missing data is printed is assumed to be represented by NaN;
%
% Table includes a column at right with annual total (for pcp), or
% annual mean (for tmp), and a row at bottom with long-term
% monthly means and long-term annual mean.  NaN's in any monthly
% value for a year yield NaN for that years total or average.  But
% the long-term monthly means and annual means are computed based
% on all available data, ignoring NaNs.
%
% fprintf uses tha \f "formfeed" character so that if more than 70
% years of data, data for years 71 and on are put on second page.

nlines=70;  % max number of years to put on page 1


% Check input mtx
[mY,nY]=size(Y);
if(nY~=13)
	error('ppt matrix should be 13 cols')
end

blnk=blanks(70);  % make 1 x 70 string of blanks

% Scale and round so that ascii file is in correct units;
% Transpose because fprintf will grab data down cols
tenfact=eval(['1e',num2str(digs)]);  % scaling factor; 
% dig=2 gives 100, the usual for TRL ppt data
D=Y(:,2:13) * tenfact;  % scale the monthly data

% Make vector of annual totals or means
% Annual value will be NaN if any monthly value NaN
if dtype==1 | dtype==4; % pcp data or pan evap data
	a=(sum(D'))';
elseif dtype==2; % tmp data
	a = ((sum(D'))')/12;
else
	error('So far, only dtype==1or2 works')
end


% Vector of long-term monthly means, based on all valid
% (non NaN) data for the month
L1=isnan(D);
na1 = sum(L1); % rv of number of NaNs in each month
na2 = sum(na1); % total number of NaNs in all months

% If any data missing, replace them with zeros, adjust the
% sample size, and compute the column means based on the
% non-NaN data only
mmy = mY(:,ones(12,1)); % constant rv of number of years in data mtx
if na2>0; % at least one missing value
	F=D; % copy of data
	F(L1)=zeros(na2,1); % replace NaNs with zeros
	ssize=mmy - na1; % sample size disregarding NaNs
	sum1 = sum(F); % column sums
	x = round(sum1 ./ ssize);  % col mean 
	% get mean of annual sums
	L2=~isnan(a);
	aa=a(L2);
	asize=sum(L2);
	xx=round(mean(aa));
else; % no NaNs, easier computation
	F=D;
	x=round(mean(F)); % long-term means
	xx=round(mean(a));
end



Z = round([Y(:,1) D a]);

% Build format for a line of printout
fmt=['%4d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%6d\n'];
fmt2=['Ave %5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%6d\n'];
fmt3=['N   %5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%6d\n\n'];
fmt4='%s\n';

% Build line 1 of header
stn=deblank(stn);
ste=deblank(ste);
id=deblank(id);
str1=[stn ',  ' ste '   ('  id ')'];
lt=length(str1); % Number of characters in string 1
if lt > 70
	error('Make title 70 or fewer characters')
elseif lt==70
	str1=str1;  % use title as is
else;  % center title
	ndiff=70-lt;
	nn=fix(ndiff/2);
	remmy = rem(ndiff,2);
	if remmy==0
		str1 = [blanks(nn)  str1 blanks(nn)];
	else
		str1=  [blanks(nn+remmy) str1 blanks(nn)];
	end
end

% Build line 2 of header
if dtype==1; % pcp
	ss1='MONTHLY PRECIPITATION - ';
elseif dtype==4; % pan evap
	ss1='MONTHLY PAN EVAPORATION - ';
else
	error('pcp or pan evap only valid data types yet');
end
ss2=deblank(dunits);
str2 = [ss1 ss2];
lt=length(str2); % Number of characters in string 1
if lt > 70
	error('Make title 70 or fewer characters')
elseif lt==70
	str2=str2;  % use title as is
else;  % center str2
	ndiff=70-lt;
	nn=fix(ndiff/2);
	remmy = rem(ndiff,2);
	if remmy==0
		str2 = [blanks(nn)  str2 blanks(nn)];
	else
		str2=  [blanks(nn+remmy) str2 blanks(nn)];
	end
end

% Build Line 3 of header
% Recall that z(1),z(2) are long and lat in decimal degrees, and
% that z(3) is elevation in m
% Want to print long-lat in deg, minutes;  elev in ft

% Convert 
londeg=round(abs(fix(z(1))));
lonmin=round((abs(z(1))-londeg) * 60);
latdeg=round(abs(fix(z(2))));
latmin=round((abs(z(2))-latdeg) * 60);
elft=round(z(3)*3.28); % meters to feet
str3a=['LATITUDE: ' int2str(latdeg) ' '  int2str(latmin) 'N     '];
str3b=['LONGITUDE: ' int2str(londeg)  ' '  int2str(lonmin) 'W     '];
str3c=['ELEVATION: '  int2str(elft) ' ft'];
str3=[str3a str3b str3c];

% Make continuation header
str4=[stn ', '  ste  '  (CONT)'];

% Make column headings 
hdra='Year   J    F    M    A    M    J    J    A    ';
hdrb='S    O    N    D Annual';
hdr = [hdra hdrb];

% Recall that the series has mY years of data.
% Compute whether 1 or 2 pages needed
if mY>nlines;
	npages=2;
else
	npages=1;
end


% Print first page

% Print lines 1-3 
fprintf(fid1,'%s',str1);
fprintf(fid1,'\n\n');
fprintf(fid1,'%s',str2);
fprintf(fid1,'\n\n');
fprintf(fid1,'%s',str3);
fprintf(fid1,'\n\n');

% print col header
fprintf(fid1,'%s',hdr);
fprintf(fid1,'\n\n');

nn= min([nlines mY]);

% Print data for up to first 70 years
for i= 1:nn;
	fprintf(fid1,fmt,Z(i,:));
end


% If only 1 page, print final two lines
if npages==1;
	fprintf(fid1,'%\n');
	fprintf(fid1,fmt4,blanks(4));
	fprintf(fid1,fmt3,ssize,asize);
	fprintf(fid1,fmt2,x,xx);
	fprintf(fid1,'\f');
else; % 2 pages, go to next page, put header info, and finish up
	fprintf(fid1,'\f');
	fprintf(fid1,'%s',str4);
	fprintf(fid1,'\n\n');
	fprintf(fid1,'%s',hdr);
	fprintf(fid1,'\n\n');
	for i=(nn+1):mY;
		fprintf(fid1,fmt,Z(i,:));
	end
	fprintf(fid1,fmt4,blanks(4));
	fprintf(fid1,fmt3,ssize,asize);
	fprintf(fid1,fmt2,x,xx);
	fprintf(fid1,'\f');
end


