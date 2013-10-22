function flwcon1
% flwcon1: convert monthly CA DWR web-site natural flow to monthly .mat file
% flwcon1:
% Last revised 12-12-00
%
% Convert monthly flow data from California Dept. of Water Resources web page
% format to .mat file monthly format
%
%
%*** INPUT 
%
% No args.
% User prompted to point to point to .txt input file
%
%
%*** OUTPUT
%
% No args.
% User prompted to point to .mat ouput file to create or overwrite
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%
%*** NOTES
%
% Infile assumed to be csv, with disregarded original (non-natural) flow between commas ignored and
% data coded as in the file obtained from CA DWR web page 11-00
%
% Output data stored as monthly tsm X in target .mat file.  Year as  col 1. Jan-Dec
% data in AF as cols 2:13
%
% Missing data written as NaN
% Leading or trailing all-missing years lopped off
% Any internal all-missing years remain in outfile
%
% Units of outfile are AF
% 

% Get infile and load data
[file1,path1]=uigetfile('*.txt','Infile of web-format monthly data');
pf1=[path1 file1];

% Read data into cell array of strings
F = textread(pf1,'%s','delimiter','\n','whitespace','');

% Find row that has "FLOW", and delete thru that
F1=char(F);
F1(:,1:6)=[]; % strip so "FLOW" is at left
i1 = strmatch('FLOW',F1);
if isempty(i1);
    error('No FLOW line in file');
end;
F(1:3)=[];
F=char(F);
[mF,nF]=size(F);
% F is now char matrix with 3 header lines deleted


% REMOVE ANY ALL-BLANK LINES FROM F
B=blanks(nF);
B=repmat(B,mF,1);
L=F==B;
i1=find((all(L'))');
if ~isempty(i1);
    F(i1,:)=[];
end;


% FIND COLUMNS OF F1 WITH ',' AND MAKE SURE SAME COLS FOR ALL ROWS
i1=findstr(F(1,:),',');
ncomma=length(i1);
if ncomma~=2;
    error('Must have 2 commas');
end;
icomma=i1;
for n = 1:ncomma;
    c=F(:,i1(n));
    L = c==',';
    if ~all(L);
        error(['Commas # ' int2str(n) ' not in same col in all rows']);
    end;
end;


% STORE YEAR, MONTH, NATURAL FLOW VECTORS
yr = str2num(F(:,1:4));
mon = str2num(F(:,5:6));
ion=icomma(2)+1;
ioff=nF;
dstr= F(:,ion:ioff);


% REPLACE ANY MISSING VALUES WITH NaN
misser = 'm  '; % will replace these with NaN
i1=strmatch(misser,dstr)
if ~isempty(i1);
    nfill = size(i1,1); 
    fillx = repmat('NaN',nfill,1);
    dstr(i1,1:3)=fillx;
end;


% DATA FROM STRING TO DOUBLE
d=str2num(dstr);


% ALLOCATE FOR DATA STORAGE
ndiff=diff(yr);
if any(ndiff<0) | any(ndiff>1);
    error('Year increments incorrectly');
end;
if any(mon<1) | any(mon>12); 
    error('Month must be between 1 and 12');
end;


% PAD START AND END IF NEEDED
yrgo = min(yr);
yrsp=max(yr);
mongo = mon(1);
monsp = mon(length(mon));
if mongo ~=1;
    mon1pad = (1:(mongo-1))';
    n1pad = length(mon1pad);
    mon=[mon1pad; mon];
    d = [repmat(NaN,n1pad,1); d];
    yr = [repmat(yrgo,n1pad,1); yr];
    
end;
if monsp ~=12;
    mon2pad = ((monsp+1):12)';
    n2pad = length(mon2pad);
    yr = [yr;  repmat(yrsp,n2pad,1)];
    mon = [mon; mon2pad];
    d = [d; repmat(NaN,n2pad,1)];
end;


% BUILD TIME SERIES MATRIX
nsize=length(d);
if mod(nsize,12)~=0; 
    error('d not divisible by 12');
end;
nrow = nsize/12;
ncol=12;
D=(reshape(d,12,nrow))';
year = (yrgo:yrsp)';


% STRIP OFF ANY ALL-NAN YEARS
i1 = find(all((isnan(D))'));
D(i1,:)=[];
year(i1)=[];


% ADD BACK ANY INTERNAL ALL-NAN YEARS
dyear = diff(year);
if any(dyear~=1);
    yrgo = min(year);
    yrsp = max(year);
    nyear = yrsp-yrgo+1;
    X = repmat(NaN,nyear,12);
    irow = year-yrgo+1; % row index
    X(irow,:)=D;
    D=X;
    year = (yrgo:yrsp)';
end;

    
% TACK ON YEAR
X=[year  D];


% Out file
[file2,path2]=uiputfile('*.mat','Outfile .mat file to hold  monthly flow data X');
pf2=[path2 file2];
eval(['save ' pf2 ' X;']);

