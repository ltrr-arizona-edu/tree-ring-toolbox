function  seaspt3(nfiles,yrs,months,k)
% seaspt3(nfiles,yrs,months,k)
%
% Seasonalize monthly climate data from several stations and store
% the seasonalized data for a common period in a matrix
%
% D Meko 9-1-95
%
%
%**************** IN ARGS *********************************
% nfiles (1 x 1)i number of monthly ppt files to get
% yrs (1 x 2)i start and end year of target matrix
% months (1 x 2)i start, end month of "season"
% k (1 x 1)i data-type: 1=ppt,2=temp (ppt is total, Temp is ave)
%
% SPECIAL: must have 13-col matrices with monthly climate data int
% to in a directory window
%
%******************** OUT ARGS ***************************
%
% No output args, but ...
%
% .MAT-FILE PRODUCED:
%
% User is prompted for a .mat filename.  In the file are stored:
%
% Y -- a time-series matrix of the seasonalized climate variable
%		for years yrs  and each station (cols)
% yr -- year cv for Y
% names -- a string matrix of input-file names corresponding to
%		the series in the columns of Y
%
%
%*********************** USER-WRITTEN FUNCTIONS CALLED 
%
% seaspt.m -- seasonalized climate matrix
%
%************** NOTES *************************
%
% If months(2)<months(1), season crosses year boundary, and will
% need one year before yrs(1) to get the seasonalized series fo
% to start in yrs(1).


a = NaN;
names=' ';

nyrs = yrs(2)-yrs(1)+1; % number of years in output matrix
yr = (yrs(1):yrs(2))';
Y=a(ones(nyrs,1),ones(nfiles,1));

% Years of monthly data uniformly needed
yrs1=yrs;
if months(2)<months(1),
	yrs1(1)=yrs(1)-1;
end
% yrs1 is start and end year of needed monthly data

 

k1=1;
n1=1; % counter for number of files read
while k1;
	clear Y1 X
	[filein, pathname] = uigetfile('*.mat', 'Select a PPT file', 500,500);
	eval(['load ',pathname,filein])

	% Check that Y1 exists.  If not, might be working with an HCN 
	% station that has X instead.
	i1 = exist('Y1');
	i2 = exist('X');
	if i1~=1 & i2~=1
		error(['Neither Y1 nor X in ' filein])
	end
	if i1~=1 & i2==1
		disp([filein 'has X instead of Y'])
		disp('press return to continue')
		Y1=X;
	end

	[m2,n2]=size(Y1);
	yyr = Y1(:,1);
	if yyr(1)>yrs1(1) | yyr(m2)<yrs1(2),
		error([filein,' does not cover needed years'])
	end
	
	% Cull needed years of monthly data
	L1 = yyr >= yrs1(1) & yyr <= yrs1(2);
	Y2 = Y1(L1,:);
	% Check that no missing data in this block
	if k==1; % precip
		L2=isnan(Y2);
		L3=Y2<0;
		if any(any(L2)) | any(any(L3))
			error('NaN or negative monthly PPT in key block')
		end
	elseif k==2; % temperature
	end


	F = seaspt(Y2,months(1),months(2),yrs1,k) 
	Y(:,n1) = F(:,2);
	names=str2mat(names,filein)

	if n1>=nfiles
		k1=0;
	end
	n1=n1+1
end

[fileout,pathname] = uiputfile('*.mat', 'Save  Seasonalized Data As');
eval(['save ' pathname fileout ' Y yr names'])


