function sitei
%
% DESCRIPTION : sitei
% Compute site indices from tree indices
%
% INPUTS  :  NONE
% OUTPUTS :  NONE
%
%
%  ZI and ZE: time series matrices with columns:
%
%  1. mean index
%  2. minimum tree index
%  3  maximum tree index
%	4. lower 95% confid limit on mean index for year
%	5. upper 95% confid limit on mean index for year
%  6. sample size (number of tree-index series averaged together)
%  7. number of cores

%
%
% F -- statistics on site chronologies: col 1 -std index,col 2 --res
%		index
%
%   1. sample size (number of years)
%   2. mean index
%   3. std dev index
%   4. mean sensitivity
%   5. first-order autocorrelation coef
%   6. 95% conf limit on (5)
%   7. second-order partial ac coef
%
%
% ACF -- lags0-5 acf and pacf
%
%   row 1: acf
%   row 2: pacf



%_________________________________________________


clear all

% Prompt for the .mat file name

% Get the .mat filename interactively; load file
[flmat,path1]=uigetfile('*.mat','Input .MAT filename ?');
pf1=[path1 flmat];
flold=flmat; % Will need this file name later
eval(['load ' pf1]);

% Check if the mat file exists in the workspace, and if output from treei.m in file
if ~exist('ET') | ~exist('IT'),
  error('This .mat file does not have ET and IT')
end

a=NaN;

% Fit history
Lwhen=exist('Fwhen','var')==1;
if ~Lwhen;
    Fwhen = cell(8,4);
end;
Fwhen{6,1}='sitei'; % function
Fwhen{6,3}=flmat; % infile
clear flmat;


% Build a mask to omit some tree indices if desired
if ~(exist('tmask')==1);
   tmask=[]; 
   tmask=logical(tmask);
else;
   tmask=logical(tmask);
end
tmask = treemask(Tnms,tmask);


% Use mask to pull subset of Tnms, ET,IT, ETn,ITn,ETyrs,ITyrs
% Tnms
trees = Tnms(tmask,:);
Eyrs = ETyrs(tmask,:);
Iyrs = ITyrs(tmask,:);
[m1,dum1]=size(Iyrs); % m1 is number of tree indices for building chron

% Size F, which will hold chron statistics
F = a(ones(7,1),ones(2,1));


% *********** First handle the standard indices

% Size Z, which will hold site chron
yrZ =[min(Iyrs(:,1)): max(Iyrs(:,2))]'; % col vector of years for Z 
mZ = length(yrZ); % number of years in site chronology
Z = a(ones(mZ,1),ones(7,1));  % initialize Z as NaN matrix


% Size A and B, which will be timeseries matrices of the tree
% indices and the number of cores for each tree 
A = a(ones(mZ,1),ones(m1,1)); 
B = a(ones(mZ,1),ones(m1,1));

% Get tree indices and put in A;  get number of cores and put in B
for j = 1:m1; % loop over trees (non-masked)
	i1 = Iyrs(j,3); % start index in IT
	i2 = i1 + (Iyrs(j,2)-Iyrs(j,1));
	% compute corresponding start, end indices of same in A,B
	k1 = 1 + Iyrs(j,1) - yrZ(1);
	k2 = k1 + (i2-i1);
	% Pull and put
	A (k1:k2,j) = IT(i1:i2);
	B(k1:k2,j) = ITn(i1:i2);
end

% Compute number of cores in each year; First fill NaN members of B
% with zero;  Then sum over rows;
B1=B;
L1 = isnan(B);
L1s = sum(sum(L1));
B1(L1)  = zeros(L1s,1);
Z(:,7) = (sum(B1'))';

% Compute number of tree indices in each year
A1=A;
L1 = ~isnan(A);
Z(:,6) = (sum(L1'))';

% Compute maximum and minimum tree-index in each year
% Algorithm assumes every valid tree index exceeds -1E7 or is less than
%  +1E7
 toobig = 1E7;
 toosmall=-1E7;
L2 = ~L1;  % logical pointer to NaNs in A
L2s = sum(sum(L2));
% First do maximum
A1=A;
A1(L2) = toosmall(ones(L2s,1),:);
Z(:,3) = (max(A1'))';
% Then minimum
A1=A;
A1(L2) = toobig(ones(L2s,1),:);
Z(:,2) = (min(A1'))';


% Mean index (chronology index) for years with  fewer than 6 trees.
% In this case, use simple arithmetic mean
% 
L3 = Z(:,6)<6 & ~isnan(Z(:,6));
C = A(L3,:);
L4 = isnan(C);
L4s = sum(sum(L4));
C(L4) = zeros(L4s,1);
zn = Z(L3,6); % number of trees
zs = (sum(C'))';
z = zs ./ zn;
Z(L3,1)=z;


% Mean index in each year and approximate 95% error bars
%
% First build vector of prob points for t dist
znmax = max(Z(:,6)); % max number trees in any year
tvect = tinv(.975,1:(znmax-1));


L5 = Z(:,6)>=6 & ~isnan(Z(:,6));
C = A(L5,:);
[mC,nC]=size(C);
z1=a(ones(mC,1),:);
z4 = z1;
z5=z1;
for i = 1:mC;
	c1 = C(i,:);
	jc = ~isnan(c1);
	y = (c1(jc))';
   
	% alternative using simple arith mean
	%ymn=mean(y);
%	stdy = std(y);
%	ny = length(y);
%	sem= stdy/sqrt(ny);
%	df = ny-1;
%	mult = tvect(round(df))* sem;
%   %mult = tinv(.975,df)*sem

	%alternative using bisq mean
	[ymn,varyh,df,w,ybar,se]=bisqmean(y);
	mult = tvect(round(df))* sqrt(varyh);


	z1(i)=ymn;
	z4(i) = ymn - mult;
	z5(i) = ymn + mult;
end
Z(L5,[1 4 5])=[z1 z4 z5];

yrZI=yrZ;
ZI=Z;  % finished standard site chronology 



% *********** Second handle the AR residual indices

% Size Z, which will hold site chron
yrZ =[min(Eyrs(:,1)): max(Eyrs(:,2))]'; % col vector of years
%   for Z 
mZ = length(yrZ); % number of years in site chronology
Z = a(ones(mZ,1),ones(7,1));  % initialize Z as NaN matrix


% Size A and B, which will be timeseries matrices of the tree
% indices and the number of cores for each tree 
A = a(ones(mZ,1),ones(m1,1)); 
B = a(ones(mZ,1),ones(m1,1));

% Get tree indices and put in A;  get number of cores and put in B
for j = 1:m1; % loop over trees (non-masked)
	i1 = Eyrs(j,3); % start index in ET
	i2 = i1 + (Eyrs(j,2)-Eyrs(j,1));
	% compute corresponding start, end indices of same in A,B
	k1 = 1 + Eyrs(j,1) - yrZ(1);
	k2 = k1 + (i2-i1);
	% Pull and put
	A (k1:k2,j) = ET(i1:i2);
	B(k1:k2,j) = ETn(i1:i2);
end

% Compute number of cores in each year; First fill NaN members
% of B
% with zero;  Then sum over rows;
B1=B;
L1 = isnan(B);
L1s = sum(sum(L1));
B1(L1)  = zeros(L1s,1);
Z(:,7) = (sum(B1'))';

% Compute number of tree indices in each year
A1=A;
L1 = ~isnan(A);
Z(:,6) = (sum(L1'))';

% Compute maximum and minimum tree-index in each year
% Algorithm assumes every valid tree index exceeds -1E7 or is
% less than
%  +1E7
 toobig = 1E7;
 toosmall=-1E7;
L2 = ~L1;  % logical pointer to NaNs in A
L2s = sum(sum(L2));
% First do maximum
A1=A;
A1(L2) = toosmall(ones(L2s,1),:);
Z(:,3) = (max(A1'))';
% Then minimum
A1=A;
A1(L2) = toobig(ones(L2s,1),:);
Z(:,2) = (min(A1'))';


% Mean index (chronology index) for years with  fewer than 6
% trees.
% In this case, use simple arithmetic mean
% 
L3 = Z(:,6)<6 & ~isnan(Z(:,6));
C = A(L3,:);
L4 = isnan(C);
L4s = sum(sum(L4));
C(L4) = zeros(L4s,1);
zn = Z(L3,6); % number of trees
zs = (sum(C'))';
z = zs ./ zn;
Z(L3,1)=z;


% Mean index in each year and approximate 95% error bars
%
% First build vector of prob points for t dist
znmax = max(Z(:,6)); % max number trees in any year
tvect = tinv(.975,1:(znmax-1));


L5 = Z(:,6)>=6 & ~isnan(Z(:,6));
C = A(L5,:);
[mC,nC]=size(C);
z1=a(ones(mC,1),:);
z4 = z1;
z5=z1;
for i = 1:mC;
	c1 = C(i,:);
	jc = ~isnan(c1);
	y = (c1(jc))';
   
	% alternative using simple arith mean
	%ymn=mean(y);
%	stdy = std(y);
%	ny = length(y);
%	sem= stdy/sqrt(ny);
%	df = ny-1;
%	mult = tvect(round(df))* sem;
%   %mult = tinv(.975,df)*sem

	%alternative using bisq mean
	[ymn,varyh,df,w,ybar,se]=bisqmean(y);
	mult = tvect(round(df))* sqrt(varyh);


	z1(i)=ymn;
	z4(i) = ymn - mult;
	z5(i) = ymn + mult;
end
Z(L5,[1 4 5])=[z1 z4 z5];

yrZE=yrZ;
ZE=Z;  % Finished residual chronology


% Update history 
ctime=clock;
ctime=num2str(ctime(4:5));
dtime=date;
Fwhen{6,2}=[dtime ', ' ctime];

% Save the vectors in a .mat file
newvars1=' tmask yrZI  ZI  trees ';
newvars2=' yrZE ZE Fwhen';
newvars=[newvars1 newvars2];
nsv=menu('Save variables?','Add new variables to original file?',...
   'Store new variables in new file?','QUIT');
if nsv==1,
    Fwhen{6,4}=Fwhen{6,3};
  eval(['save ' pf1 newvars ' -append']);
elseif nsv==2,
  [ofmat,path2]=uiputfile('*.mat','New .MAT file for storing new variables: ');
  Fwhen{6,4}=ofmat; 
  pf2=[path2 ofmat];
  eval(['save ' pf2 newvars]);
end
% End of file
