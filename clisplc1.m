function clisplc1
%
% Splice station climate series
%
% Meko 1-14-97
%
%********************INFILES ***********************
%
% Two .mat files of data for the same station.  Assumed continuous
% years and NaN for missing values.
%
%****************** OUTFILE ********************
%
% A .mat file of the altered first input file.  This altered file
% now might include updated or substituted data from the second
% file.
%
%********************* NOTES *******************************
%
% Example.  Have two ppt .mat files for Bisbee.  One is the
% long-term Bradley series, which ends in 1963.  Another is
% a NWS/COOP Bisbee that starts in 1961 and ends in 1985.
%
% Want to whenever possible replace any missing values from the
% first series with data from the second.  Also want to extend
% the first data with data from the second series if the second
% series has more recent data.

a=NaN;

% Get first file -- this file is to be altered and/or extended
[file1,path1]=uigetfile('*.mat','input .mat file with base series');
pf1=[path1 file1];
eval(['load ' pf1]);
Z1=Z;
clear Z;

% Get second file -- this file contains data to be substituted in
% the first file.  Also might extend the first file
[file2,path2]=uigetfile('*.mat','input .mat file with more data');
pf2=[path2 file2];
eval(['load ' pf2]);
Z2=Z;
clear Z;

% Size
[m1,n1]=size(Z1);
[m2,n2]=size(Z2);

% Get year vector for the two input files
yr1=Z1(:,1);
yr2=Z2(:,1);

% Specify sup-period of master input series to use and
% get the data for the subperiod
k1=input('Accept entire period of input master series? [Y]/N ','s');
if isempty(k1) |  k1=='y' | k1=='Y'
	k1='Y';
elseif k1=='N' | k1=='n'
	k1='N';
	yrgo=input('Start year of subperiod of master to be used: ');
	yrsp=input('End year of subperiod of master to be used: ');
	L1=yr1>=yrgo & yr1<=yrsp;
	Z1=Z1(L1,:);
	yr1=Z1(:,1);
	[m1,n1]=size(Z1);
else
	error('Invalid input for k1')
end


% Specify subperiod of slap-on series to use and get the data
k1=input('Accept entire period of input slap-on series? [Y]/N ','s');
if isempty(k1) |  k1=='y' | k1=='Y'
	k1='Y';
elseif k1=='N' | k1=='n'
	k1='N';
	yrgo=input('Start year of slap-on data to use: ');
	yrsp=input('End year of slap-on  data to use: ');
	L2=yr2>=yrgo & yr2<=yrsp;
	Z2=Z2(L2,:);
	[m2,n2]=size(Z2);
	yr2=Z2(:,1);
else
	error('Invalid input for k1 on slap-on series')
end


% Any replacement possible?  Is if there is any overlap in
% the selected parts of Z1 and Z2
if min(yr2)<=max(yr1);
	yr3=(min(yr2):max(yr1))';
	L3=yr1>=yr3(1) & yr1<=yr3(length(yr3)); % points to master
	L4=yr2>=yr3(1) & yr2<=yr3(length(yr3)); % points to slap-on
	W1=Z1(L3,2:13);
	W2=Z2(L4,2:13);
	L3a=isnan(W1) & ~isnan(W2); % points to months with NaN in master and not slap on
   sum1=sum(sum(L3a));
   if sum1>0;
      W1(L3a)=W2(L3a);
   end;
   Z1(L3,2:13)=W2; % replace any overlap data of master with that of slap-on
end;
Z=Z1;


% Any appending to the end of the master series?
if max(yr2)>max(yr1);
	L6=yr2>max(yr1);
	Z3=Z2(L6,:);
	[m3,n3]=size(Z3);
	yrZ3=Z3(:,1);
	yr4=((max(yr1)+1):max(yr2))';
	nyrs4=length(yr4);
	% Fill series with Nan
	A=a(ones(nyrs4,1),ones(12,1));
	L7=yr4>=yrZ3(1) & yr4<=yrZ3(m3);
	A=[yr4 A];
	A(L7,:)=Z3;
	Z=[Z1; A];
end


[file3,path3]=uiputfile('*.mat','Revised master .mat file');
pf3=[path3 file3];
pf3=strtok(pf3,'.');
eval(['save ' pf3 ' Z']);


fclose all
