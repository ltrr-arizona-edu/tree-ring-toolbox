function respfun0
% respfun0:  pre-processing function to build input for respfun1.m
% respfun0;
% Last revised: 7-11-99
%
% Builds a .mat storage file of required input for response function analysis via 
% respfun1.m.  User is prompted for various types of information and is prompted to
% click on files (e.g., file with tree-ring index, file with monthly precipitation data).<P>
%
% Refer to the input instructions to respfun1.m for a definition of inputs that are 
% prepared by respfun0.  
%
%*** INPUT
% 
% No input args
% User prompted for lots of information:
% -click on a file with tree-ring index
% -whether that index is a standard, residual, or arstan chronology
% -a one-letter labelling code for first type of climate variable (e. g., P)
% -click on the file with the first climate variable
% -a one-letter labelling code for second climate variale (e.g., T)
% -click on the file with the second climate variable
% -the ending month for the plausible monthly climate window of influence (e.g., October)
% -number of months in the plausible monthly climate window (e.g., 14)
% -number of seasons to combine monthly climate data into
% -season label, start month, and end month for each season 
%  ( e.g.,   Summer    July September)
% -order of AR model for optional whitening of tree-ring series (e.g., 2)
%    (this model is fit to the climate-calibration period, and is used to
%     ensure that the tree-ring variable is not autocorrelated for that period)
% -whether or not to apply this model to AR filter the tree-ring series before
%   proceeding with correlations and response functions
% -period for analysis (e.g., 1901 1990)
% -whether climate variables are to be summed or averaged over months in seasonalizing
% -whether plots are to be color or B/Wa
% -the name of the .mat file the input information for respfun1.m is to be stored in
%
% The string matrix strout stored in your output .mat file summarizes all the 
% settings you have specified.
%
%*** OUTPUT
%
% A storage file with information needed for running respfun1.m
% 
%*** REFERENCES -- none
%*** UW FUNCTIONS CALLED
%
% nonano1
% several sub-functions
%
%*** TOOLBOXES NEEDED -- none
%
%*** NOTES
%
% Some menus have not-yet-implemented items



close all; % close any open figure windows


%*****************  SET UP THE COLOR TRIPLETS FOR COLOR AND B/W PLOTS
%
% This section hard-coded.  Normally no need to change it.  Computer nerds
% may want to play with changes. The triplets control the density of shading
% for confidence bands.  Rows 1-3 give colors for 95% 90% and other for
% "P" variable;  rows 4-6 for "T" variable;  Rows 7-9 give B/W colors that
% apply to both P and T if B/W option used

Clr = [0 0 1; .7 .7 1; 1 1 1; 1 0 0; 1 .7 .7; 1 1 1; .2 .2 .2; .6 .6 .6; 1 1 1];
set1 = ['Clr '];

%************  TREE-RING SERIES

[X,strouta]=respfn0a;

%***************************************


%-------  Some user info on required climate input files
h1=helpdlg({'Ready to specify climate variables and files.',...
      'If data in ascii files, assumed space separated.',...
      'Climate data should be in 13-col matrix, with year as col 1.',...
      'First climate variable is generally P, second T'},'Notice');


%************ MONTHLY CLIMATE DATA FOR CLIMATE VARIABLE #1 (THE "P" VARIABLE)
[P,pf2,vlabelb,stroutb]=respfn0b;


%************ MONTHLY CLIMATE DATA FOR CLIMATE VARIABLE #2 (THE "T" VARIABLE)
[T,pf3,vlabelc,stroutc]=respfn0c;

%------- combine the labels
vlab={vlabelb,vlabelc}; %Labels for climate variables (e.g., P and T)
vlab=char(vlab);

%************************  MONTHS FOR FULL MONTHLY ANALYSIS    
monnames = {'January','February','March','April','May','June','July',...
      'August','September','October','November','December'};
% Make cell of first 3 letters of each month
monnms3=char(monnames);
monnms3=monnms3(:,1:3);
monnms3=cellstr(monnms3);
[begmo,endmo,nmos,yrcross,stroute]=respfn0e(monnames,monnms3);

%**************  SEASONS
[nseas,slab,I1,stroutd]=respfn0d(monnames,begmo,endmo,nmos,yrcross); 


%************  MAXIMUM POSSIBLE AR ORDER FOR PREWHITENING
[prewh,qmax,stroutf]=respfn0f;


%************** COMPUTE START AND END OF FEASIBLE YEAR COVERAGE FOR CORRELATION
%
% Look at time coverage of X, P, T, and find longest continuous period with
% no NaNs.  For P ant T, this would be no NaNs in any months of the year

yrX = X(:,1); yrP = P(:,1); yrT = T(:,1);
yrsx=[min(yrX) max(yrX)]; yrsp =[min(yrP) max(yrP)]; yrst=[min(yrT) max(yrT)];

% Make a year vector that covers full lengths of all series
yrmin = min([yrX; yrP; yrT]);
yrmax = max([yrX; yrP; yrT]);
yrall = (yrmin:yrmax)';
nall = length(yrall);

% Make same-length logical vectors, initialized to zero, for X, P, T
LX = zeros(nall,1);  LP = zeros(nall,1); LT=zeros(nall,1);

% Make pointers to full coverage of X , P, T in yral
LX1 = yrall>=yrsx(1) & yrall<=yrsx(2);
LP1 = yrall>=yrsp(1) & yrall<=yrsp(2);
LT1 = yrall>=yrst(1) & yrall<=yrst(2);

% Make logical vectors telling if data in each year has no NaN
LX2 =  (all((~isnan(X))'))'; % cv==1 if no NaN in the year
if strcmp(prewh,'Yes') & qmax>0; % might need leading years of tree ring data
   itemp=find(LX2);
   itempgo = itemp(1);
   LX2(itempgo:(itempgo+qmax-1)) = 0;
end
LP2 =  (all((~isnan(P))'))'; % 
LT2 =  (all((~isnan(T))'))'; % 
if strcmp(yrcross,'Yes');
   itemp=find(LP2);
   itempgo=itemp(1);
   LP2(itempgo)=0;
   itemp=find(LT2);
   itempgo=itemp(1);
   LT2(itempgo)=0;
end


% Substitute
LX(LX1)=LX2;
LP(LP1)=LP2;
LT(LT1)=LT2;
LX=logical(LX); LP=logical(LP); LT=logical(LT);

% Mark all years with valid data in all three series
Lall = LX & LP & LT;

% Compute the longest continous stretch of valid data overlapping the 3 variables
Lvalid=nonan01(Lall);
yrvalid = yrall(Lvalid);
yrsval = [min(yrvalid) max(yrvalid)]
stroutc1='MAXIMUM POSSIBLE MODELING PERIOD FOR GIVEN DATA, SEASONS AND AR ORDER';
stroutc1=char(stroutc1,...
   ['  ' int2str(yrsval(1)) '-' int2str(yrsval(2))]);


%**************  INTERACTIVE CHOICE OF ANALYSIS PERIOD
[yrs,stroutc2]=respfn0h(yrsval);


%************ MISCELLANEOUS OPTIONS
[kopt,stroutg,vlab]=respfn0g(vlab);

%************** BUILD CALL TO RESPFUN1.M
call='respfun1(P,T,X,yrs,kopt,nmos,endmo,prewh,qmax,vlab,I1,slab,Clr);';


%*************  COMBINE TEXT INFO
strout=char(strouta,stroutf);
strout=char(strout,stroutb);
strout=char(strout,stroutc);
strout=char(strout,stroutc1);
strout=char(strout,stroutc2);
strout=char(strout,stroute);
strout=char(strout,stroutd);
strout=char(strout,stroutg);


%*****************  SAVE THE WORKS
set1=' X P T kopt nmos endmo prewh qmax vlab slab I1 yrs Clr call strout '; 

[file6,path6]=uiputfile('*.mat','Store respfun1.m inputs in this mat file');
pf6=[path6 file6];
eval(['save ' pf6 set1 ';']);







%************ SUBFUNCTIONS ********

%******** TREE RINGS

function [X,strout]=respfn0a
%respfn0a:  subfunction of respfun0 that gets the tree-ring chronology

kmen1 = menu('Type of Source File for Tree-Ring Series',...
   'ITRDB-style .crn', ...
   '2-column ascii, with year in col 1',...
   'Matlab chrononology file',...
   '.mat file with 2-column matrix');
kchron = menu('Chronology type',...
   'Standard',...
   'Residual',...
   'ARSTAN');
switch kchron;
case 1;
   chrontyp='Standard Chronology';
case 2;
   chrontyp ='Residual Chronology';
case 3;
   chrontyp = 'ARSTAN Chronology';
end

%-------- Get and store tree-ring data
switch kmen1;
case 1; % .crn file
   [file1,path1]=uigetfile('*.crn','Input .crn file with tree-ring chron');
   pf1 = [path1 file1];
   [x,s,yr]=crn2vec2(pf1); % call user-written function to get the data
   % Lop-off pre-1850 data because not needed for response function, generally
   Ltemp = yr<1850;
   if any(Ltemp);
      x(Ltemp)=[]; yr(Ltemp)=[];
   end
   % Store the data
   X=[yr x];
   % Cleanup
   clear x s yr Ltemp;
case 2; % 2-column ascii file
   [file1,path1]=uigetfile('*.dat','Input .dat file with tree-ring chron');
   pf1 = [path1 file1];
   eval(['X=load(pf1);']);
   Ltemp = X(:,1)<1850;
   if any(Ltemp);
      X(Ltemp,:)=[];
   end
case 3;  % chronology .mat storage file
   [file1,path1]=uigetfile('*.mat','Input .mat file with tree-ring chron');
   pf1 = [path1 file1];
   % Chronologt will be first col in ZI or ZE, depending on if standard or
   % residual chron, and year will be in yrZI or yrZE
   if kchron==1; % standard chron
      eval(['load ' pf1 ' ZI yrZI;']);
      yr = yrZI;
      x = ZI(:,1);
      clear ZI yrZI;
   elseif kchron==2; % residual chron
      eval(['load ' pf1 ' ZE yrZE;']);
      yr = yrZE;
      x = ZE(:,1);
      clear ZE yrZE;
   else;
      error('kchron must be 1 or 2 for chrons in .mat storage files');
   end
   Ltemp = yr<1850;
   if any(Ltemp);
      x(Ltemp)=[]; yr(Ltemp)=[];
   end
   X = [yr x];
   
case 4; % .mat file with year in column 1, data in col 2
   [file1,path1]=uigetfile('*.mat','Input .mat file with tree-ring chron');
   pf1 = [path1 file1];
   % Chronology will be in matrix X, with year as col 1
   eval(['load ' pf1 ' X;']);
   yr = X(:,1);
   x = X(:,2);
   clear X;
   Ltemp = yr<1850;
   if any(Ltemp);
      x(Ltemp)=[]; yr(Ltemp)=[];
   end
   X = [yr x];
end

strout='TREE-RING SERIES';
strout = char(strout,['  File = ' pf1]);
strout=char(strout,['  Version = ' chrontyp]);


%**************** FIRST CLIMATE VARIABLE

function [P,pf2,vlabel,strout]=respfn0b
% respfn0b: subfunction for respfun0.m.  First monthly climate data.

% 1-letter code for type of climate variable
prompt={'Enter 1-letter code for first climate variable'};
def = {'P'};
titdlg = 'First Climate Variable';
lineNo=1;
answer=inputdlg(prompt,titdlg,lineNo,def);

vlabel=answer{1};

kmen1 = menu(['Choose file-type for input ' vlabel ' climate variable'],...
   'Ascii file with 13-col matrix',...
   '.mat file holding 13-col matrix');
switch kmen1;
case 1; % P data in ascii file
   [file2,path2]=uigetfile('*.dat',['Input .dat file with ' vlabel ' climate data']);
   pf2 = [path2 file2];
   eval(['P=load(pf2);']);
   clear file2 path2;
case 2; % P data in .mat file
   [file2,path2]=uigetfile('*.mat',['Input .mat file with ' vlabel ' climate data']);
   pf2 = [path2 file2];
   eval(['S = load(pf2);']);
   F = fieldnames(S);
   if ~isempty(strmatch('X',F,'exact'));
      P=S.X;
   elseif ~isempty(strmatch('Y',F,'exact'));
      P=S.Y;
   elseif ~isempty(strmatch('Z',F,'exact'));
      P=S.Z;
   else
      error([pf2 ' does not have X, Y or Z']);
   end
   clear file2 path2 F S;
end

strout=['CLIMATE VARIABLE: ' vlabel];
strout=char(strout,['  File of monthly data: ' pf2]);


%*************** SECOND CLIMATE VARIABLE


function [T,pf3,vlabel,strout]=respfn0c
% respfn0c: subfunction for respfun0.m.  Second monthly climate data.

% 1-letter code for type of climate variable
prompt={'Enter 1-letter code for second climate variable'};
def = {'T'};
titdlg = 'Second Climate Variable';
lineNo=1;
answer=inputdlg(prompt,titdlg,lineNo,def);

vlabel = char(answer{1});

kmen1 = menu(['Choose file-type for input ' vlabel ' climate variable'],...
   'Ascii file with 13-col matrix',...
   '.mat file holding 13-col matrix');
switch kmen1;
case 1; % T data in ascii file
   [file2,path2]=uigetfile('*.dat',['Input .dat file ' vlabel ' climate data']);
   pf3 = [path2 file2];
   eval(['T=load(pf3);']);
   clear file2 path2;
case 2; % T data in .mat file
   [file2,path2]=uigetfile('*.mat',['Input .mat file with ' vlabel ' climate data']);
   pf3 = [path2 file2];
   eval(['S = load(pf3);']);
   F = fieldnames(S);
   if ~isempty(strmatch('X',F,'exact'));
      T=S.X;
   elseif ~isempty(strmatch('Y',F,'exact'));
      T=S.Y;
   elseif ~isempty(strmatch('Z',F,'exact'));
      T=S.Z;
   else
      error([pf3 ' does not have X, Y or Z']);
   end
   clear file2 path2 F S;
end

strout=['CLIMATE VARIABLE: ' vlabel];
strout=char(strout,['  File of monthly data: ' pf3]);


%**************  SEASONS


function [nseas,slab,I1,strout]=respfn0d(monnames,begmo,endmo,nmos,yrcross) 
% respfn0d: subfunction for respfn0.m; specify seasons

prompt={'Enter number of seasons'}; def = {'4'};
titdlg = 'Number of seasons for seasonal analyis';
lineNo=1;
answer=inputdlg(prompt,titdlg,lineNo,def);
nseas = str2num(answer{1});
strout='SEASONS';
strout = char(strout,['Number of seasons = ' int2str(nseas)]);

% Build info text
txt1a = ['Full Window (' int2str(nmos) ' mos) is '];
txt1b = [char(monnames(begmo)) '-' char(monnames(endmo))];
txt1 = [txt1a txt1b];

I1 = repmat(NaN,nseas,2); % to store index of start and end months for seasons
% these are relative to the nmos in the full monthly window

% Build a cell variable with month names for all nmos 
if strcmp(yrcross,'No'); % if not ross year boundary
   monnms = monnames(begmo:endmo);
else;
   monnms=[monnames(begmo:12) monnames(1:endmo)];
   nprev = 12-begmo+1; % number of months in year t-1
end

%  Build cell variable with 3-letter starts of month names
monnms3=char(monnms);
monnms3=monnms3(:,1:3);
monnms3=cellstr(monnms3);

lastsp=0;
icross=0; % year boundary not yet crossed
nrun = 0 ; % running total of months in seasons specified so far

% Loop over seasons
for n = 1:nseas;
   
   % Set relative index for first month of season
   if n==1;
      mongo = 1;
   else;
      mongo = lastsp+1;
   end
   
   % User specify name and ending month of season
   mongonm = monnms{mongo}; % name of starting month
   monspnm = ''; % name of stop month
   prompt={'Season code','Start Month','End month'};
   def = {['Seas' int2str(n)],mongonm,monspnm};
   lineNo= 1;
   titdlg = ['Season #' int2str(n) ' of ' int2str(nseas) '; ' txt1];
   answer=inputdlg(prompt,titdlg,lineNo,def);
   slab{n}=answer{1}; % store season name   
   
   % Compute the month "number" for ending month
   monspnm = answer{3};
   istop = strmatch(upper(monspnm(1:3)),upper(monnms3));
   if isempty(istop);
      error([monspnm ' is invalid stop month']);
   else;
      if n==2;
         n;
      end;
      
      % If two matches for stop month, ask whether in year t or t-1
      if  length(istop)>1;
         strstop=menu({['which ' monspnm ' do you mean for end month '],...
               [' of ' slab{n}]},...
               [monspnm ' in year t-1'],...
               [monspnm ' in year t']);
         if strstop==1; 
            istop=istop(1);
         else;
            istop=istop(2);
         end;
      end;
        
      
      monspnm=monnms{istop};
      if length(istop)>1;
         if icross==0 & istop(1)>lastsp; % have not yet crossed year boundary
            lastsp = istop(1);
            if lastsp>=nprev; % you reached the year boundary with this seas
               icross=1;
            else;
            end
         else
            lastsp = istop(2);
            icross=1;
         end
      else;
         lastsp = istop;
      end
   end
   
   I1(n,:)= [mongo lastsp];
   
   stradd1 = int2str(n);
   stradd2 = char(slab{n});
   stradd3 = [mongonm '-' monspnm];
   stradd = ['  ' stradd1 ' ' stradd2 ' = ' stradd3];
   
   strout=char(strout,stradd);      
end

if lastsp~=nmos;
   error('Seasons do not cover all nmos months in month window');
end

icross;

%******** MONTHS

function [begmo,endmo,nmos,yrcross,strout]=respfn0e(monnames,monnms3)
% respfn0e: subfunction for respfun0.m.  Months


prompt={'Enter ending month for climate window of monthly analysis'};
def = {'September'};
titdlg = 'Ending Month';
lineNo=1;
answer=inputdlg(prompt,titdlg,lineNo,def);
endname = answer{1};

% Check that first three chars of endname match those of a month
imonth = strmatch(upper(endname(1:3)),upper(monnms3),'exact');

if isempty(imonth);
   error([endname ' is invalid month']);
else
   endmo = imonth;  % index of ending month
   endname=monnames{imonth}; % set to fully-spelled month name
end


prompt={'Enter number of months in window for monthly analysis'};
def = {'14'};
titdlg = 'Number of Months';
lineNo=1;
answer=inputdlg(prompt,titdlg,lineNo,def);
nmos = str2num(answer{1});

% Compute whether monthly window crosses year boundary
if nmos>endmo;
   yrcross='Yes';
else
   yrcross='No';
end


%------- Compute beginnning month for full monthly window
if strcmp(yrcross,'Yes');
   begmo = endmo-nmos+13;
else
   
   begmo=endmo-nmos+1;
end
begname = monnames{begmo};


strout=['MONTHS'];
strout=char(strout,['  Number of months in full window = ' int2str(nmos)]);
strout=char(strout,['  Starting Month = ' begname]);
strout=char(strout,['  Ending Month = ' endname]);


%************* AR FILTERING OF SITE CHRON

function [prewh,qmax,strout]=respfn0f
% respfn0f: subfunction of respfun0.m.  Set max allowable AR model for 
% AR filtering tree-ring series before climate analysis

prewh= questdlg('Filter tree-rings by the best-fit AR model before correlation?');
if strcmp(prewh,'Cancel');
   prewh='No';
end
if strcmp(prewh,'Yes');
   txt1='  Site chronology AR-filted before correlating with climate';
elseif strcmp(prewh,'No');
   txt1='  Site chronology not AR-filtered before corrrelating with climate';
end

if strcmp(prewh,'Yes');
   prompt={'Enter AR order'};
   def = {'3'};
   titdlg = 'Maximum-order AR model to consider in filtering calibration-period tree rings';
   lineNo=1;
   answer=inputdlg(prompt,titdlg,lineNo,def);
   qmax = str2num(answer{1});
else
   qmax=0; % "zero-order" model means no prewhitening
end;


strout = 'AR MODELING OF CHRONOLOGY FOR CLIMATE-OVERLAP PERIOD';
strout = char(strout,txt1);
if strcmp(prewh,'Yes');
   strout=char(strout,['  Highest-order candidate AR model = AR(' int2str(qmax) ')']);
end


%********* MISC OPTIONS


function [kopt,strout,vlab]=respfn0g(vlab)
% respfn0g:  subfunction of respfun0.m.  Set options.
%

strout='MISCELLANEOUS OPTIONS';

kopt=repmat(NaN,1,4);

%-------   sum or average over months to get seasonalized climate variables

kmen =menu('How to compute seasonalized climate from monthly',...
   'Sum over months for P, average for T',...
   'Sum over months for both P and T',...
   'Average over months for both P and T');
kopt(1)=kmen;
strout=char(strout,...
   '  How to seasonalize climate variables');
if kmen==1;
   strout=char(strout,...
      '    Sum over months for P, average for T');
elseif kmen==2;
   strout=char(strout,...
      '    Sum over months for both P and T');
elseif kmen==3;
   strout=char(strout,...
      '    Average over months for both P and T');
end

   
%----------  method for confidence limits on correlation function

%kmen=menu('Method for confidence limits on correlation function',...
 %  'Theoretical -- null hypoth of zero-correlated, normally distrib populations',...
 %  'Bootstrap -- tree rings resampled',...
 %  'Circular -- tree-ring series made circular and shifted');
%if kmen~=1;
 %  error('Only option allowed so far is theoretical');
%end

%strout=char(strout,...
%   '  Method for confidence limits on correlation function');
%if kmen==1;
%   strout=char(strout,...
 %     '    theoretical -- Panofsky');
%elseif kmen==2;
%   strout=char(strout,...
 %     '    Bootstrap -- tree rings resampled');
%elseif kmen==3;
%   strout=char(strout,...
 %     'Circular -- tree-ring series made circular and shifted');
%end
%kopt(2)=kmen;
kopt(2)=1;


%-------  Option for color vs B/W plots

kmen =menu('Color or B/W Plots?',...
   'Color',...
   'B/W');
kopt(3)=kmen;
strout=char(strout,...
   '  Color or B/W Plots');
if kmen==1;
   strout=char(strout,...
      '    Color');
elseif kmen==2;
   strout=char(strout,...
      '    B/W');
end


%---------- Option for removal of cross-correlation between the 2 climate variables

%kmen = menu('Climate variable uncoupling choice',...
 %  'None: climate variables P and T',...
 %  'Climate variables P, T.P',...
 %  'Climate variables P.T, T');

%strout=char(strout,...
 %  '  Climate variable uncoupling?');
%if kmen==1;
%   strout=char(strout,...
%      '    none');
%   vlab=vlab;
%elseif kmen==2;
 %  strout=char(strout,...
  %    '    P, T.P');
  % strtemp = [vlab(2,:) '.' vlab(1,:)];
  % vlab = char(vlab(1,:),strtemp);
%else
   %strout=char(strout,...
    %  '    P.T, T');
   %strtemp = [vlab(1,:) '.' vlab(2,:)];
   %vlab = char(strtemp,vlab(2,:));
%end
%kopt(4)=kmen;
kopt(4)=1;




%*************** ANALYSIS PERIOD

function [yrs,strout]=respfn0h(yrsval)
% respfn0h:  subfunction of respfun1.m .  years for analysis


yrs=repmat(NaN,1,2);

txt1='Period for analysis (must be within ';
txt2=[int2str(yrsval(1)) '-' int2str(yrsval(2)) ')'];
txt3=[txt1 txt2];
txt4=['Start year earlier than ' int2str(yrsval(1)) ' or end year later than ' int2str(yrsval(2))];


kwh1=1;

while kwh1==1;
   prompt={'Enter start year','Enter ending year'};
   def = {int2str(yrsval(1)),int2str(yrsval(2))};
   titdlg = txt3;
   lineNo=1;
   answer=inputdlg(prompt,titdlg,lineNo,def);
   yrs(1)=str2num(answer{1});
   yrs(2)=str2num(answer{2});
   
   kwh1=0;
   
   if yrs(2)<=yrs(1);
      h1=helpdlg('Start year must precede end year','Hit OK and try again!');
      pause(3);
      if ishandle(h1);
         close(h1);
      end
      kwh1=1;
   end
   
   if yrs(1)<yrsval(1) | yrs(2)>yrsval(2);
      h1=helpdlg(txt4,'Hit OK and try that again');
      pause(3);
      if ishandle(h1);
         close(h1);
      end
      kwh1=1;
   end
   
end

strout='SELECTED PERIOD FOR ANALYSIS OF CLIMATE-TREE RING RELATIONSHIP';
strout=char(strout,...
   ['  ' int2str(yrs(1)) '-' int2str(yrs(2))]);   





