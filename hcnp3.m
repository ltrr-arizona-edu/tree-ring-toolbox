function ns=hcnp3
% ns=hcnp3
%
% D Meko 5-6-96
%
% Replace missing-values in "filnet" version of HCN data 
% pnet2.dat with "original" data from porg2.dat if non-missing
% data exists for that station/year/month in porg2.dat
%
%****************** INPUT ***********************************
%
% No arguments.  hcnp3 needs data from several files, whose
% file names are hard coded:
%
% porg2.dat -- 'original' version of hcn data
% pinvorg.dat -- corresponding station info 
% pnet2.dat -- 'filnet' version of hcn data
% pinvnet.dat -- corresponding station info
%
% The above 4 files were produced by two runs of hcnp2.for
%
%
%************************ OUTPUT ******************************
%
% hcnp3.m has only one output argument:
%
% ns (1 x 1)i  -- number of station
%
% The main output of hcnp3.m is the "revised" version of the filnet
% monthly data.  This revised version is stored in the matrix X.
% The user is prompted for a .mat filename to save X in.
%
%
%********************* NOTES ************************************
%
% hcnp3.m runs fast after the initial slow loading of the
% huge .dat files porg2.dat and pnet2.dat.  Loading takes about
% 1'10'', but the whole function runs in about 2'.



a=NaN;
xmiss=-99.99;

disp('Starting to load porg2.dat and pnet2.dat')
if(~exist('porg2'))
	eval(['load porg2.dat'])
end
if(~exist('pnet2'))
	eval(['load pnet2.dat'])
end
disp('Finished loading those files')


[mX,nX]=size(pnet2);

disp('Initializing storage for X')
X = a(ones(mX,1),ones(nX,1));

fid3=fopen('pinvnet.dat','rt');
fid4=fopen('pinvorg.dat','rt');


% Read all lines in file3 to get number of stations
ns=0;
while 1
	line = fgetl(fid3);
	ns=ns+1;
 	if ~isstr(line), break, end
% 	disp(line)
end
ns=ns-1;  % number of stations
disp(['Total of ns = ',int2str(ns),' stations'])

frewind(fid3);


% Loop over stations
for n=1:ns;
	disp(['Starting on station ',int2str(n)])
	lnf=fgetl(fid3); % filnet
	lno=fgetl(fid4); % original
	idf=lnf(6:11);
	ido=lno(6:11);
	% Get start,end years and row indices, filnet
	yrfgo=str2num(lnf(38:41));
	yrfsp=str2num(lnf(43:46));
	rowfgo=str2num(lnf(50:54));
	rowfsp=str2num(lnf(58:62));
	% Get start,end years and row indices, original
	yrogo=str2num(lno(38:41));
	yrosp=str2num(lno(43:46));
	rowogo=str2num(lno(50:54));
	rowosp=str2num(lno(58:62));

	% make a NaN matrix same row-size as the filnet series 
	nyrsf = yrfsp-yrfgo+1;
	yrf = (yrfgo:yrfsp)';
	A = a(ones(nyrsf,1),ones(14,1));

	% Pull the row-subset of original data
	O=porg2((rowogo:rowosp),:); 
	yro=O(:,2);	

	% Make row-index vector to slip original data into
	% correct rows of A
	i1 = yro-yrf(1)+1;
	A(i1,:)=O;	

	% Pull rows of filnet series
	F=pnet2((rowfgo:rowfsp),:);
	Fold=F;
	F2=F(:,3:14);

	A2=A(:,3:14);

	% Logical pointer to missing values that can be replaced
	L1= F2==xmiss;
	L2= A2~=xmiss;
	L3=L1&L2;

	% Replace missing values
	ss=sum(sum(L3));
	F2(L3)=A2(L3);
	F(:,3:14)=F2;

	% Convert all missing values to NaN
	L1 = F2==xmiss;
	if any(any(L1));
		ss = sum(sum(L1));
		F2(L1) = a(ones(ss,1),:);
	end

	% Put month columns of filnet back with id and year cols
	F(:,3:14)=F2;

	% Store result for this series in X
	X(rowfgo:rowfsp,:)=F;

end

disp('hcnp.mat is recommended name for the next .mat file')
files=uiputfile('*.mat','Save the changed filnet F in this file');
eval(['save ',files,' X'])

fclose('all')
