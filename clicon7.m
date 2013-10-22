function Z=clicon7(pf1)
% Convert an NCDC Coop network monthly pcp file to matlab file
% Z=clicon1(fnm)
% Last revised 1-19-00
% 
% NCDC Coop monthly pcp conversion.  Source is an ascii file obtained by 
% downloading from NCDC. Desired form is 13-col matrix Z, with year as 
% first column, other cols as monthly pcp in inches, missing values as 
% NaN, and no internal rows missing.
%|
%*** IN
%
% pf1 -- optional path\filename of file with the input monthly pcp as
%   obtained from NCDC. If not supplied, you are prompted to click on file.
%
%
%*** OUT   
%
% Z (mZ x nZ)r   monthly pcp matrix, with year as col 1
%
% You are also prompted for a .mat file to save this matrix Z in
%
%*** REFERENCES -- none
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED  -- NONE

xmiss=-99999; % missing value code in inpu

if nargin ~=1;
   [file1,path1]=uigetfile('*.txt','Infile of Coop NCDC data');
   pf1=[path1 file1];
else;
end;


% Read coop file in as a text array
A=textread(pf1,...
   '%s','delimiter','\n','whitespace','');
A=char(A);

% Pull year column
year = str2num(A(:,18:21));
nyr = length(year);

% Compute input columns for monthly data
igo=35:12:167; % starting col for monthly data
isp=igo+5;

% Initialize storage matrix
X=repmat(NaN,nyr,13);
X(:,1)=year;

% Get the data and store it
for i=1:12;
   i1=igo(i);
   i2=isp(i);
   c = str2num(A(:,i1:i2));
   X(:,i+1)=c;
end;


% Replace any missing values
L=X==xmiss;
if any(any(L));
   X(L)=NaN;
end

% Check that years in order
if any(diff(year)<=0);
   error('Years out of order or duplicate');
end;

% Check that no entire missing years; if there are, make em NaN
if any(diff(year)~=1);
   yrgo=min(year);
   yrsp=max(year);
   nnyr = yrsp-yrgo+1;
   yrnew = (yrgo:yrsp)';
   Z = repmat(NaN,nnyr,13);
   
   islot = year-yrgo+1; % row in Z to put each year of data from X
   Z(:,1)=yrnew;
   Z(islot,:)=X;
else;
   Z=X;
end;

% Scale data from hundredths of inches to inches
Z(:,2:13)=Z(:,2:13)/100;


% Output file
[file2,path2]=uiputfile('*.mat','Output mat file to store Z in');
pf2=[path2 file2];
eval(['save ' pf2 ' Z;']);
   
   
   
