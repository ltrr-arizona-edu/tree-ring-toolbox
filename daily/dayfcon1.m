function Y=dayfcon1
% 
% Convert USGS "netscape" file of daily flow in cfs to a 
% monthly matrix of acre-ft, NaN filled
%
% Meko, 8-17-96
%
%
%****************** IN ********************************
%
% Assumes 
% - you pulled in a monthly data file over the internet
% - you used edix to cut off the header info and convert the
% -  '.' to blank so that file is all numeric
% - file now has 4 col
%		year
%		month
%		day
%		discharge (cu-ft/sec)



% Get the data file
[file1,path1]=uigetfile('*.dat','Daily netscape flow file to work on');
pf1=[path1 file1];
eval(['load ' pf1]);


% Put file contents in matrix D
f1=strtok(file1,'.');
eval(['D = ' f1 ';']);


% Store vectors of years, months, days, and cfs
yr=D(:,1);
mon=D(:,2);
day=D(:,3);
d=D(:,4);


% Find min and max years, and size storage matrix Y
yr1=min(yr);
yr2=max(yr);
nyrs=yr2-yr1+1;
a=NaN;
Y = a(ones(nyrs,1),ones(13,1));
Y(:,1)=(yr1:yr2)';

% Compute the conversion factor.  This times cfs is acre-ft/day
f=(24*3600)/43560;


% Loop over years. Check that some data exists for each year.
% 
for n=yr1:yr2;
	irow=n-yr1+1; % row index for this year in the output matrix
	L1=yr==n;
	s1=sum(L1);
	if s1==0,
		disp(['No data for year ' int2str(n)]);
		error('Abandon ship!')
	end
	disp(['Working on year ' int2str(n)]);
	% Loop over months.  Check that some data exists for each month.
	for j=1:12;
		L2= L1 & mon==j;
		s2=sum(L2);
		if s2==0,
		 disp(['No data for year/month ' int2str(n) '/' int2str(j)]);
		end	
		d1=d(L2); % the daily data
		% Sum the daily values for this year/month and convert
		Y(irow,j+1)=f * sum(d1);
	end
end


% SAVE THE MONTHLY ACRE-FT TOTALS IN A .MAT MATRIX
[file2,path2]=uiputfile('*.mat','Save monthly acre-ft mtrix Y here');
pf2=[path2 file2];
eval(['save ' pf2 ' Y']);
