function clisplc3
% clisplc3: splices recent monthly climate data on older existing file
% clisplc3;
% Last revised 2-19-00
%
% Used to update files of monthly precipitation or temperature, while giving 
% information on agreement of the recent and older data for overlapping years.
%
%*** INPUT -- no args
%
% You are prompted to point to input .mat files with the old and recent (update) data.
% These files should have data in 13-column monthly climate format. The year is 
% column, Jan-Dec data is columns 2-13.  Units should be the same in the two files.
% The files should not be missing any internal years.  In other words, rows should 
% be continuous years.
%
% You are prompted for:
% Type of data -- P, T, or PDSI
% Whether to overwite existing data values with new data (NOTES).
% Whether to accept overwriting an old data value with
%   a recent-data NaN (NOTES)
%
%
%*** OUTPUT -- no args
%
% Graphical output includes 2 fiure windows. These summarize the closeness of the update 
% data and the old data for overlapping years (NOTES).
%
% You are prompted for the file to write the updated data to. This can be the 
% original old data file, or a new file. 
%
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% Graphical Output.  
% Window 1 plots differences between the monthly data in old and recent files along 
% with box plots showing distribution of monthly values in the old file. This lets 
% you judge the diffences in the context of the variability of the old data, individually 
% year by year, for the most recent up to 4 years of overlap.
% Window 2 gives a text table comparing the maximum difference of old and recent 
% data with the standard deviation and range of the old data.
%
% Overwriting existing dat with new data.  Because the recent data overlaps the old data, 
% there may be some months/years that overwriting old with recent will change existing 
% old data.  You might not want to do this.  For example, if the old data is HCN adjusted 
% and the new data is unadjusted COOP data, you may want to keep your old data up to the 
% point when you need data to extend the file to present.  To illustrate:
% 
% Old file might be HCN 1900-1995
% Recent file is COOP, 1993-1999, with 1993-95 included to allow overlap for comparison
% Blanket overwriting of old with recent would replace the HCN 1993-95 data with the
%  COOP data, which is probably a bad idea
%
% But be carefule with this. Sometimes the last values of the old file might have been  
% provisional at the time, and the more recent data 'better.'  It is best to look at the
% updating station by station.  You can also manually pare the new data file down to 
% restrict the data that could possibly be overwritten. 
%
% NaN handling.  Sometimes the recent data may have NaNs in months/years that 
% the old file has data for.  This might happen, for example, if the old file had 
% had data estimation done, or if the old file had been constructed on a different rule 
% of how many daily values could be missing before marking the monthly data as missing.
% You are prompted for how to handle overwriting of 'good' data with NaN:
%    1 - Do it.  It is OK to overwrite a data value with a NaN.
%    2 - Do not do it.  Do not overwrite any existing old data values with NaNs. 
% 

close all;

% PROMPT FOR DATA TYPE
kmen2=menu('Choose Date Type',...
   'Precipitation',...
   'Temperature',...
   'PDSI');
switch kmen2;
case 1; % P
   units='(in)';
   dtype='P';
case 2; % T
   units='(\circF)';
   dtype='T';
case 3; %PDSI
   units='';
   dtype='PDSI';
end;
txt1='\Delta';
ylab =[txt1 dtype units];


% GET FILE NAMES AND INPUT DATA

[file1,path1]=uigetfile('*.mat','Input old data file -- to be updated');
pf1=[path1 file1];
[file2,path2]=uigetfile('*.mat','Input recent-data file -- to update with');
pf2=[path1 file2];

% Get old data, and store in X
eval(['s=load(''' pf1 ''');']) % old data, in either X,Y or Z
f=fieldnames(s);
Lx=strcmp(f,'X');
Ly=strcmp(f,'Y');
Lz=strcmp(f,'Z');
Lall =[Lx; Ly; Lz];
if ~any(Lall);
   error([pf1 ' does not contain X Y or Z']);
end;
if sum(Lall)>1;
   error([pf1 ' contains more than one of X Y or Z']);
end;
if Lx;
   X=s.X;
elseif Ly;
   X=s.Y;
elseif Lz;
   X=s.Z;
else;
   error('X not assigned from structure s');
end;
clear ans s Lx Ly Lz f Lall;

% Get recent data, and store in Y
eval(['s=load(''' pf2 ''');']) % old data, in either X,Y or Z
f=fieldnames(s);
Lx=strcmp(f,'X');
Ly=strcmp(f,'Y');
Lz=strcmp(f,'Z');
Lall =[Lx; Ly; Lz];
if ~any(Lall);
   error([pf2 ' does not contain X Y or Z']);
end;
if sum(Lall)>1;
   error([pf2 ' contains more than one of X Y or Z']);
end;
if Lx;
   Y=s.X;
elseif Ly;
   Y=s.Y;
elseif Lz;
   Y=s.Z;
else;
   error('Y not assigned from structure s');
end;
clear ans s Lx Ly Lz f Lall;

txttit=['Updating ' file1 ': Difference New-Old in Overlap'];
   


% FIND OVERLAP YEARS.  RECALL X IS OLD DATA, Y RECENT.

% Get year vectors for X and Y
[mX,nX]=size(X);  [mY,nY]=size(Y);
yrX=X(:,1); yrY=Y(:,1);
yrgoX=yrX(1);  yrspX=yrX(mX);
yrgoY=yrY(1);  yrspY=yrY(mY);
if yrspX>yrspY;
   error('Old data ending year is later than New data ending year');
end;


% Mark rows of X overlapping Y
LX=yrX>=yrgoY & yrX<=yrspY;
X1 = X(LX,:); yrX1=X1(:,1);  % the subset of old data

% Mark rows of Y overlapping X, w/o regard to NaNs
LY= yrY>=min(yrX1) & yrY<=max(yrX1);
Y1=Y(LY,:);  yrY1 = Y1(:,1); % the subset of recent data

% Mark years with dual data in at least one month
X1a = X1(:,2:13); Y1a=Y1(:,2:13);
L1=~isnan(X1a) & ~isnan(Y1a);
L2=(any(L1'))'; % cv marking years with some common data
if ~any(L2);
   error('No months/years with data in common');
end;
clear L1;

% Compute and mark years with comparison data.  Only consider most recent up-to-4 years
nsum2=sum(L2);
if nsum2<=4;
   L3=L2;
else;
   i2=find(L2);
   i2=i2((nsum2-3):nsum2);
   L3=repmat(0,length(L2),1);
   L3=logical(L3);
   L3(i2)=1;
   clear i2;
end;
years1 = yrX1(L3);  % overlap years to go in figure legend
nyears1 = length(years1);

% Cull overlap data, with NaNs considered
X2=X1(L3,:);  Y2=Y1(L3,:);
X2a=X2(:,2:13); Y2a=Y2(:,2:13);

% Logical pointer to comparison data for each month
L4=~isnan(X2a) & ~isnan(Y2a);  % col 1 is for Jan, ...


% COMPUTE MONTHLY DIFFERENCES RECENT MINUS OLD, AND STORE INFO IN CELLS
dep=cell(1,12); % will hold differences
D=repmat(NaN,12,nyears1);

%rv of # of comparison values in each month of year
if size(L4,1)==1;
   nsum4=L4;
else;
   nsum4=sum(L4);  
end;

for n =1:12;
   LL4=L4(:,n);
   nn=nsum4(n);
   if nn==0; % no comparison values for this month
      dep{n}=NaN;
      D(n,:)=NaN
   else;
      dep{n}=Y2a(LL4,n) - X2a(LL4,n);
      D(n,LL4)=Y2a(LL4,n)' - X2a(LL4,n)';
   end;
end;

% RECALL
% dep{} is row cell of recent-old values
% years1 is cv of years
% cols of L4 mark which years the values in dep{} apply to


% COMPUTE STATS FOR OLD DATA
Xa=X(:,2:13);
L6=~isnan(Xa);
nsum6=  sum(L6);
if any(nsum6<2);
   error('Must have at least 2 years of old data');
end;
s=nanstd(Xa);  rng = range(Xa); % standard dev and range
xamed=nanmedian(Xa); % monthly medians
xam  = repmat(xamed,size(Xa,1),1);


% BOXPLOTS  OF OLD DATA, AS DEPARTURES FROM MONTHLY MEDIANS
figure(1);
boxplot(Xa-xam);
ylabel(ylab);
hold on;

% SUPERPOSE COMPARISON DATA
sym={'o','*','d','v'};
for n = 1:nyears1;
   d=D(:,n); % cv
   plot((1:12)',d,sym{n});
end;

hold off;
title(txttit);
grid;
zoom xon;

%build legend
stryr=num2str(years1);
ylimm=get(gca,'YLim');
ydel=diff(ylimm)/20;
ypt=[ylimm(2)-ydel; ylimm(2)-2*ydel; ylimm(2)-3*ydel;  ylimm(2)-4*ydel];
xpt=9.3;
for n=1:nyears1;
   text(xpt,ypt(n),[sym{n} ' '  stryr(n,:)]);
end;


% FIGURE 2 -- TEXT INFO
figure(2);
str1=sprintf('%7.0f',nsum6);
if size(D,2)==1;
   dmax=D;
else;
   dmax=(max(D'))';
   dmin=(min(D'))';
   Lneg=abs(dmin)>abs(dmax);
   if any(Lneg);
      dmax(Lneg)=dmin(Lneg);
   end;
end;
j=(1:12);
D = [j; nsum6; dmax'; s; rng]
S=sprintf('%3.0f %3.0f %6.2f %6.2f %8.2f\n',D);
title('Maximum departure (new-old) in context of variability of monthly data');
stra='M  N   dmax    s     rng';
text(0.1,.8,stra);
text(0.1,0.75,S,'VerticalAlignment','top');
set(gca,'visible','off');

% BUILD STORAGE MATRIX FOR UPDATED FILE
% Series will begin in yrgoX and end in yrspY
nyr = yrspY-yrgoX+1;
XX = repmat(NaN,nyr,13);
yrXX = (yrgoX:yrspY)';
XX(:,1)=yrXX;
% Put old data in
Lx = yrXX>=yrgoX & yrXX<=yrspX;
XX(Lx,:)=X;


% PROMPT FOR USER OPTION -- WHAT YEAR TO START REPLACING OLD WITH NEW DATA
kmen3=menu('Choose One',...
   'Start updating in first year after end of old data',...
   'Start replacing with first year of available new data',...
   'Specify year to start updating');
if kmen3==1;
   yron = yrspX+1;
elseif kmen3==2;
   yron = yrgoY;
elseif kmen3==3;
   prompt={'Enter the first year:'};
   def={num2str(yrspX+1)};
   dlgTitle='First year to start replacing or updating';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);
   yron  = str2num(answer{1});
   if yron>(yrspX+1);
      error ('first year for updating cannot be later than first yr after end of X');
   end;
   if yron<yrgoY;
      error('First year for updating cannot be earlier than first year of new data');
   end;
end;

yroff = yrspY; % last year for updating is last year of new data

% PULL BLOCK OF DATA FROM XX TO BE REPLACED
Lxx = yrXX>=yron & yrXX<=yroff;
U =XX(Lxx,:);
Ua = U(:,2:13);
Lgood = ~isnan(Ua); % mark non-NaN data in block to be replaced

% PULL BLOCK OF DATA FROM NEW FILE -- this will  be the replacing data
Lyy = yrY>=yron & yrY<=yroff;
V = Y(Lyy,:);
Va = V(:,2:13);

% PROMPT FOR USER OPTION -- NAN OVERWRITING
kmen1=menu('Choose NaN Overwriting Option',...
   'Overwrite only NaNs in the old file',...
   'Overwrite existing data as well as NaNs in the old file');
if kmen1==1;
   Va(Lgood)=Ua(Lgood); % replace the new data in replacement block with old data
   V(:,2:13)=Va;
elseif kmen1==2;
end;

% REPLACE & UPDATE & rename
XX(Lxx,:) = V;
Z = XX;


% WRITE OUTPUT
kmen4=menu('Choose one',...
   ['Overwrite ' file1 ' with revised data'],...
   'Write revised data to new file',...
   'Bail out ... need to think this over');
if kmen4==1;
   eval(['save ' pf1 ' Z -append']);
elseif kmen4==2;
   [file3,path3]=uiputfile('*.mat','Output file to store revised data');
   pf3=[path3 file3];
   eval(['save ' pf3 ' Z ;']);
elseif kmen4==3; 
   % no action
end;


% 
