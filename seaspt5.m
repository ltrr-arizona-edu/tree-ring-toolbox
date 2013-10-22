function seaspt5
%
% Seasonalize multi-station matrix of regcli2.m output climate
% data.  Make matrix for one season. Each row a year, each col
% seasonalized data for a station.  Written to process regcli2.m
% output.
%
% Meko 2-8-98
%
%********************** INPUT FILES  & PROMPTS **********************
%
% User points to a data file previously produced by regcli2.m. This
% file has year in col 1, monthly data for region 1 in cols 2-13,
% for region 2 in cols 14-25, etc.
%
% User also prompted for desired:
% 
% begmo: start month of season (1=jan, 12=dec)
% endmo:  end month...
% kopt: 1=ppt, 2=tmp  (note: ppt is summed, tmp averaged over months)
%
%********************* OUTPUT FILE *****************************
%
% User points to file to store output. This is a .mat file with
% year in col 1 and seasonalized data for regions in the other
% columns.  Both a .mat and .dat version are produced


% Open the input file
[file1,path1]=uigetfile('*.mat','Input monthly data file');
pf1=[path1 file1];
eval(['load ' pf1]);


% ppt or tmp
kopt=input('Precip (1) or Temp (2)?  ');
% Months to start, end season
begmo=input('Start month of season: (1,...12):  ');
endmo=input('End month of season:  (1,...,12):  ');

% Compute number of regions in the matrix.  Easy because 12 columns
% per region assumed, with leading year col. Data assumed in X.
yr = X(:,1); % store year vector
yrs=[min(yr) max(yr)]; % start, end years of monthly data
X(:,1)=[]; % remove year column from X
[mX,nX]=size(X); % 
nreg=nX/12;  % number of regions

% Allocate for multi-region mtx to store seasonalized data
% Trial call to ptseas to get correct row size for seasonalized mtx
W=[yr X];
F=seaspt(W,begmo,endmo,yrs,kopt);
yrseas=F(:,1);
[mF,nF]=size(F);
a=NaN;
Y=a(ones(mF,1),ones(nreg,1));
clear F;

% Loop over regions
for n=1:nreg;
	jgo=1+ (n-1)*12; % Jan data this station in source mtx
	jsp=jgo+11;  % Dec data ...
	jcols=jgo:jsp;
	W=[yr X(:,jcols)]; % 13-col matrix, monthly data this region
	% Seasonalize
	F=seaspt(W,begmo,endmo,yrs,kopt); % year in col 1, seas data in col 2
	% Store data for this region
	Y(:,n)=F(:,2);
end

% Remove all NaN years if needed
L1= isnan(Y);
L2=  (all(L1'))';
if sum(L2)>0;
	Y(L2,:)=[];
	yrseas(L2)=[];
end

% add the year col
Y=[yrseas Y];
nrows=size(Y,1);

% .mat file of seasonalized data
[file2,path2]=uiputfile('*.mat','File for output seasonalized mtx');
pf2=[path2 file2];
eval(['save ' pf2  ' Y']);

% .dat file of same
blnk=' ';
flout=[strtok(pf2,'.')  '.dat'];
fid2=fopen(flout,'w');
fmt1='%4.0f';
fmt2='%6.2f';
fmt3='%s\n';
for n = 1:nrows;
	year=Y(n,1);
	y=Y(n,2:(nreg+1));
	fprintf(fid2,fmt1,year);
	fprintf(fid2,fmt2,y);
	fprintf(fid2,fmt3,blnk');
end
fclose (fid2);



