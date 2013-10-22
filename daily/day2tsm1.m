function [X,f,stn] = day2tsm1(S);
% day2tsm1: block of daily flow data, USGS style, to a time series matrix
% [X,f,stn] = day2tsm1(S);
% Last revised 2-15-02
%
% Utility function called by dayusgs1.m to reorganize daily flow data into a time series matric
%
%
%*** INPUT
%
% S (? x ?)s   character matrix of USGS-formated daily data (see dayusgs1.m)
%
%
%*** OUTPUT
%
% X (? x 5)r   daily time series matrix: year, month, day-of-month, Julian day, and daily flow value in cfs
% f (? x 1)s   streamflow-value qualification code (see USGS Web site)
% stn (1 x ?)s   USGS staion number (see notes)
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED -- NONE
%
%*** NOTES
%
% stn.  A check is made that the same USGS station number is in each line of S.  If not, an error is returned


[mS,nS]=size(S);


s=S(1,:);
sfirst=s;
slast=S(mS,:);
if ~strcmp(s(1:4),'USGS');
    error('chars 1-4 of first line not USGS');
end;
if ~strcmp(slast(1:4),'USGS');
    error('chars 1-4 of last line not USGS');
end;


i1 = find(isspace(s));
i3=i1-1;
i2=[1 i1+1];

i2(end)=[];
S1=S(:,i2(1):i3(1)); % USGS

% USGS station number
S2=S(:,i2(2):i3(2)); 
% Check that the id for all lines is the same
s2=S2(1,:);
SS2=repmat(s2,mS,1);
Lcheck = S2==SS2;
if ~all(all(Lcheck));
    error(['Some row of S has a different station number than ' s2]);
end;
stn = s2;
    

Sdate = S(:,i2(3):i3(3)); 

disp('Converting year to numeric');
Syr = Sdate(:,1:4);
yr =str2num(Syr);
disp('Converting month to numeric');
Smonth = Sdate(:,6:7); 
monthe= str2num(Smonth);
disp('Converting day to numeric');
Sday = Sdate(:,9:10); 
day=str2num(Sday);
N=size(day,1);

% Compute first and last year for target daily tsm; 
yrgo = yr(1);
yrsp =yr(end);
nyrX = yrsp-yrgo+1;

% Allocate tsm
nlines = nyrX*366; % number of lines
X=repmat(NaN,nlines,5);

%----  BUILD FIRST 4 COLS OF X

% Build calendar year vector
i1 = 0:1:(nyrX-1);
I1=repmat(i1,366,1);
J1= yrgo + I1;
yrvec = J1(:); % year vector
clear i1 I1 J1 ;
% Build month vector;
j=[repmat(1,31,1); repmat(2,29,1);  repmat(3,31,1);  repmat(4,30,1); ...
        repmat(5,31,1);  repmat(6,30,1);  repmat(7,31,1);  repmat(8,31,1);...
        repmat(9,30,1);  repmat(10,31,1);  repmat(11,30,1); repmat(12,31,1)];
J=repmat(j,1,nyrX);
mon = J(:);
clear J;
% Build day of month vector
j=[1:31   1:29   1:31 1:30 1:31 1:30 1:31 1:31 1:30 1:31 1:30 1:31]';
J=repmat(j,1,nyrX);
daymon= J(:);
clear J;
% Build day vector
j=(1:366)';
J=repmat(j,1,nyrX);
dayseq = J(:);

X(:,1:4)=[yrvec mon daymon dayseq];


disp('Processing and storing the vectors of daily values and flags -- hang on awhile');
S(:,1:i3(3))=[]; % now has only the data and the flag
[mS,nS]=size(S);
x=repmat(NaN,mS,1); % to store vector of data
f = repmat(blanks(1),mS,1); % to store flag for estimation
jday = repmat(NaN,mS,1); % to store julian day
for n =1:mS;
    d=S(n,:); % a daily flow , and possibly an estimation flag
    
    L1 = ~isspace(d);
    dd=diff(L1);
    ilast=find(dd==-1);
    ifirst = find(L1);
    ifirst=ifirst(1);
    ione = sum(dd==1); % 1 if just a value, 2 if also a flag
    if ione==2;
        ddd = [0 dd];
        i2 = find(ddd==1);
        ff = d(i2(2));
        f(n)=ff;
    else;
    end;
    clear ddd i2 
    x(n)=str2num(d(ifirst:ilast));
    jday(n)  = mdy2jdy(monthe(n),day(n));
end;

disp('Inserting the daily data into X, the time series matrix');

% Recall:  yr, monthe and day are the year, month and day of data in x

% Compute insertion index

islot = ((yr- yrgo)*366)  +  jday;
X(islot,5)=x;

A=[];