function wy2mat01
% wy2mat1:  monthly Calif DWR water year precip spreadsheet to .mat file
% wy2mat1;
% Last revised 3/27/01
%
% Used to convert a slew of .mat precip files derived from .xls files from Maury Roos 
% into my standard .mat monthly climate files (see notes)
% 
%*** INPUT
%
% No args
% Function treats all .mat files in the current directory.  Assumes those files
% have W, elev, lat, lon, id, name.
% User prompted for output directory, which should be empty.  
%
%
%*** OUTPUT
%
% New .mat files, with monthly data in tsm X, and structure S:
%  S.lat, S.lon, S.elev, S.id, S.name
% Overall summary file roosp1.mat, with
%   D.lat, D.lon, D.elev, D.filename, D.id, D.name, D.first, D.last
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% Maury sent me ~ 74 .xls station ppt files, as .xls spreadsheets.  I used excel and copy/paste
% to put the data from these into individual .mat files, with 
% W (14 cols) holding 1) water year, 2) water year sum, 3-14) Oct-Sept values, in inches
% elev in feet
% lat, lon, in decimal deg
% id: identification code
% name: station name
%
% Problems. 
% W has 14 columns, not 13.  Do not want sum as col.  
% W is organized by water year,
%   with Oct, t-1, as col 3, Sept, year t as col 14.  I want Jan-Dec.  
% W may have gaps of years without data.  I want the year col to be continous
%
%
%  STEPS
%  PROMPT FOR SOME SETTINGS
%  BUILD CELL STRING OF THE INPUT .MAT FILENAMES; ALLOCATE


%  PROMPT FOR SOME SETTINGS

% Output directory
prompt={'Enter directory:'};
def={'c:\projs\aj1\data\nsierrap\calendar\'};
dlgTitle='Output directory';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
path2=answer{1};


%  BUILD CELL STRING OF THE INPUT .MAT FILENAMES; ALLOCATE

D= dir('*.mat'); % structure.  D.name has file name
nfile=size(D,1); % number of files
Data.lat=repmat(NaN,nfile,1); % latitude
Data.lon=repmat(NaN,nfile,1);
Data.elev=repmat(NaN,nfile,1);
Data.filename=cell(nfile,1);
Data.stnname=cell(nfile,1);
Data.first=repmat(NaN,nfile,1);
Data.last=repmat(NaN,nfile,1);
Data.units='inches';
Data.how='wy2mat01.m';


for n = 1:nfile; % loop over stations
    dfn = D(n,1).name; % input file name
    disp(['Starting on ' dfn]);
    pf2=[path2 dfn]; % output path/file 
    eval(['load ' dfn ';']);  % load data for station
    % Assign location and identifying info
    Data.lat(n)=lat;
    Data.lon(n)=lon;
    Data.elev(n)=elev;
    Data.filename{n}=dfn;
    Data.stnname{n}=name;
    
    yrW=W(:,1);
    [mW,nW]=size(W);
    
    
    % Flesh out W to continuous years
    yrgo = W(1,1);
    yrsp = W(mW,1);
    yrV = (yrgo:yrsp)';
    nyrV = length(yrV);
    V=repmat(NaN,nyrV,14);
    V(:,1)=yrV;
    irow = yrW-yrgo+1;
    V(irow,:)=W;
    
    % Trim all except monthly cols from V
    V = V(:,3:14);
    
    % Reorganize V into calendar year as X
    B = V(:,1:3); % oct-dec
    A = V(:,4:12); % jan-sept
    
    yrX =[(yrV(1)-1); yrV];
    A1=[(repmat(NaN,1,9)); A];
    B1=[B; (repmat(NaN,1,3))];
    X=[A1 B1];
    [mX,nX]=size(X);
    if all(isnan(X(1,:)));
        yrX(1)=[];
        X(1,:)=[];
        [mX,nX]=size(X);
    end;
    if all(isnan(X(mX,:)));
        yrX(mX)=[];
        X(mX,:)=[];
        [mX,nX]=size(X);
    end;
    X=[yrX X];
    Data.first(n)=yrX(1);
    Data.last(n)=yrX(mX);
    
    clear S;
    S.name = name;
    S.id = id;
    S.lat=lat;
    S.lon=lon;
    S.elev=elev;
    S.first=yrX(1);
    S.last=yrX(mX);
    S.units='in';
    vwhat='Calendar-year organized precip for N Sierra Nevada, from CA DWR';
    vwhat=char(vwhat,'Manually manipulated from .xls files, then converted via wy2mat01.m');
    S.what = vwhat;
    eval(['save ' pf2 ' X S;']);
    disp(['    Finished saving ' pf2]);
        
    
end;
eval(['save ' path2 'Xinfo Data;']); % save all station id and location info in file Xinfo.mat
    
   
disp('here')
 
    