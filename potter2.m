function potter2(Z,yr,P,alph)
%
% Potter statistic (see potter.m) applied to a matrix of time series
%
% Meko 8-4-96
%
%*********************  IN ARGS **************************************
%
% Z (mZ x nZ)r tsm, no year col;  mZ variables, nZ years;
%		any col may have leading or trailing NaN
% yr (mZ x 1)i   years vector for Z
% P (mZ x ?)i pointer matrix telling which cols of Z are to be
%		averaged to form the "true" climate series for testing the
%		individual station series
%		Each row is for a time series of Z.  Thus, for example,
%
%		P =[2 4 7 9 0 0 0 0;
%          1 6 5 0 0 0 0 0...
%
%		says that series 1 is to be tested agains the 4-station
%		mean series of stations 2,4, 7 and 9; and that
%		station 2 against  the mean series for stations 1,6, and 5, etc
%
% alph (1 x 1)r  significance level for T to be used in tabular
%		summary of passes and fails
%
%
%************************ IN FILES *********************************
%
% User is prompted for name of an ascii file with labels for 
% stations in Z.  
%
%************************ FILE OUTPUT ********************************
%
% User is prompted for name of a .txt file to hold information on the 
% station groupings and a table of results. Contents are like this
%
%
%                     Potter Summary
%
% Station Identification -- key stations and their grouped true clim
%
% 1 Benson		7 8 9 
% 2 Bisbee		1 4 8 12
%
%
% Test results
%
%   STATION	MAX T	YEAR	SIGNIF	SHIFT(IN)	FLAG
%
% 1 Benson		4.9	1950	.01	+0.40    	**
% etc
%
%
% 
%***********************************************************

global tablepot

yrgo=min(yr);
yrsp=max(yr);
yrs1=yrsp-yrgo+1;
YRS1=[yrgo yrsp;yrgo yrsp];

% Get filename with station identifier input info
[file1,path1]=uigetfile('*.txt','Input station filenames');
pf1=[path1 file1];
fid1=fopen(pf1,'r');

% Get the filename for output ascii info
[file2,path2]=uiputfile('*.txt','Output ascii info here');
pf2=[path2 file2];
fid2=fopen(pf2,'w');



% Allocate

[mZ,nZ]=size(Z); % matrix of nZ time series


for n=1:nZ; % loop over stations
	v=Z(:,n); % questionable  series
	
	% Get station file or name tag for questionable series
	% Put in a 10-character string
	cf=fgetl(fid1);
	cf = cf(~isspace(cf));
	ncf=length(cf);
	cfb=blanks(10);
	if ncf>10;
		cf=cf(1:10);
		ncf=10;
	end
	cfb(1:ncf)=cf;
	str7=sprintf('%10s ',cfb);


	% Form "true" series
	p=P(n,:); % pointer, zero filled at right
	n1=sum(p>0); % number of series to be averaged to form true climate
	p=p(1:n1); % truncate to get rid of trailing zeros
	W=Z(:,p);
	I1=ones(1,n1)'; % use all series in subset in call to smav1
	k1=3; % NaN option for tsmav1
	[w,yrw,nn,s]=tsmav1(W,YRS1,k1,I1);

	% Get common segment of test and true series
	[x,y,yr2]=overlap(v,w,YRS1);
	YRS2=[min(yr2)  max(yr2)];
	[T,I0,T0,dstari0,sig,T05]=potter(x,y,tablepot);

	str1=sprintf('%3.0f ',n);
	str2=sprintf('%5.1f ',T0);
	str3=sprintf('%5.1f ',T05);
	yr0=yr2(I0); % shift after this year
	str5=sprintf('%4.0f ',yr0);
	str4=sprintf('%5.2f ',dstari0);
	str6=sprintf('%4.0f %4.0f ',YRS2(1,:));
	strall=[str1 str7 str6 str2 str3 str5 str4];
	fprintf(fid2,'%s\n',strall);
end


fclose(fid2);
