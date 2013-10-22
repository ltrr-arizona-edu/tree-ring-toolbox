function day2fle1
% day2fle1:  Eischeid daily P file to individual daily matrix files
% dayfle1; 
% Last revised  12-4-00
%
%*** IN
%
% User prompted for pf1,path and filename of data file as obtained from John Eischeid (see notes)
% User prompted for mode of run
%   Pass1 = scan the file and record line numbers of headers and data
%   Pass2 = use the info file from pass 1 to access data and store in individual station files
%
%*** OUT
%
% No output arguments
% One of two possibilities
%   Pass1:  .txt file of station info and line numbers; also used later as access database (see notes) 
%   Pass2:  .mat files P??????.mat for individual stations  ( see notes)
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED 
% dayorg1 -- organize a month of daily obs
% leapyr
%
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% pf1:  File was got by ftp from Eischeid public space as .gz by M Munro. MM then used perl to extract
%   subsets for running this.  Note that original complete files massive -- 600 MB+.  These files have one
%   station after another.  Header line before each.  Data then:  line 1 for year 1 month 1, all days; 
%   line 2 for year 1 month 2, all days, etc.
% Pass1 and Pass2.  Day2fle1 is run twice. First to survey.  Second to store data in .mat files
% .mat files;  a station's daily data are in Y: Columns are   
%   1-year (calendar)
%   2-month (jan=1)
%   3-day of month (use 29 for Feb)
%   4-day number (1 to 366 within calendar year)
%   5-daily pcp (inches)
%   6 - code 
%
% NaN is placed in bogus days (e.g., Feb 29 of nonleap years)
% Eischeid data is purportedly "filled in" and so has no other NaNs


[file1,path1]=uigetfile('c:\data\eischeid\dlypcp\pdly?.dat','Infile of Eischeid data');
pf1=[path1 file1];


kmen1=menu('Choose 1',...
   'Pass 1: survey file and store .txt info file',...
   'Pass 2, option a: store daily data in station .mat files',...
   'Pass 2, option b: make daily time series matrix');
if kmen1==1;
   kmode='Pass1';
else;  % Note that with pass 2, option b is the only option implemented so far; thus kmen1==3 works but not kmen1==2
    if kmen1==2;
        error('I decided not to implement this option');
        return;
    end;
   kmode='Pass2';
   fid=[];
end;
fid1=fopen(pf1,'r');

t0=clock;


switch kmode;
case 'Pass1'; % survey
   ihead=repmat(NaN,2000,1);
   kwh1=1;
   nline=0;
   nfile=0;
   C=[];
   while kwh1;
      c=fgetl(fid1);
      nline=nline+1;
      if any(isletter(c));
         c(length(c))=[];
         nfile=nfile+1;
         ihead(nfile)=nline; % store line number of header line
         % Build filename
         f1=c(1:2);
         f2=c(3:6);
         k1=isspace(f1);
         if any(k1);
            f1(k1)='0';
         end;
         k1=isspace(f2);
         if any(k1);
            f2(k1)='0';
         end;
         fn=['P' f1 f2];
         c=[c ' ' fn '.mat'];
         C=char(C,c);

         disp(c);

         
      end;
      
      if feof(fid1);
         
         kwh1=0;
      end;
      
   end; % while kwh1
   
   C(1,:)=[];
   
   % Compute number of data lines (excluding header) for each station
   ihead=ihead(1:nfile);
   i1=[ihead; (nline+1)];
   num1=diff(i1)-1;
   
   % Close in file and save results
   fclose(fid1);
   [file2,path2]=uiputfile('c:\data\eischeid\dlypcp\psurvey.mat','Intermediate output to be read next pass');
   pf2=[path2 file2];
   eval(['save ' pf2 ' C  ihead num1;']);
  
case 'Pass2';  % Store daily data
   
   % Load info file stored in pass 1
   
   [file3,path3]=uigetfile('c:\data\eischeid\dlypcp\psurvey.mat','Intermediate info file produced by pass1');
   pf3=[path3 file3];
   eval(['load ' pf3 ' ihead num1 C;']);
   
   
   nfile=size(ihead,1); % number of stations in Psurvey
   
   
   if kmen1==2; % if making individual station mat files
       
       for n =1:nfile; % Loop over stations
           c=C(n,:); % char info this stn
           cnm = c(9:30); % most of stn name
           cdata=fgetl(fid1); % get line from data file
           if ~strcmp(cnm,cdata(9:30));
               disp(cnm);
               disp(cdata(9:30));
               error('Miss-match on above cnm & cdata line');
           end;
           
           %    Allocate space for daily vector
           
           mline=num1(n); % number of lines of data to read
           igo = ihead(n)+1; % first data line
           isp = igo + mline-1; % last ine
           X=repmat(NaN,31*mline,6); % more than enough space
           
           %    Build output filename
           
           ic1=findstr(c,'.mat');
           c1=c(1:ic1);
           ic=find(isspace(c1));
           filen =c((max(ic)+1):(ic1+3));
           
           % Loop over data lines
           
           ion=1; % initialize target line counter
           for m = 1:mline; % loop over lines for this station
               if m==1; 
                   kopt=1;
                   yrold=[];
                   mold=[];
               else;
                   kopt=2;
               end;
               
               cline=fgetl(fid1); % get line from data file
               [y,ycode,yrnew,mnew]=dayorg1(cline,kopt,yrold,mold);
               
               nvals = length(y);
               
               vyear = repmat(yrnew,nvals,1);
               vmon = repmat(mnew,nvals,1);
               vday = (1:nvals)';
               
               ioff=ion+nvals-1;
               X(ion:ioff,[1 2 3  5 6])=[vyear vmon vday y ycode];
               
               ion=ioff+1;
               yrold=yrnew;
               mold=mnew;
               
               
               
           end; % for m=1:mline
           
           
           %--- Trim trailing NaN rows off X
           
           L=  (all((isnan(X))'))';
           X(L,:)=[];
           [mX,nX]=size(X);
           
           
           %--- Add sequential day index (1:366)
           
           if rem(mX,366)~=0;
               error('Row size of X not evenly divisible by 366');
           else;
               n366=mX/366;
               v366=repmat((1:366)',n366,1); 
               X(:,4)=v366;
           end;
           
           % Compute information for X
           I=[X(1,[1 2 3]); X(mX,[1 2 3])];
           vlist = ['X daily eischeid time series matrix for ' cnm];
           vlist=char(vlist,'  col 1= year');
           vlist=char(vlist,'  col 2= month');
           vlist=char(vlist,'  col 3= day of month');
           vlist=char(vlist,'  col 4= 1:366 sequential day');
           vlist=char(vlist,'  col 5= daily precip in hundredths of inches');
           vlist=char(vlist,'  col 6= eischeid estimation code (see notes)');
           vlist=char(vlist,'I  row 1: first year, month, day');
           vlist=char(vlist,'   row 2: last year ,month, day');
           
           
           % Save in .mat file
           eval(['save ' filen ' X I vlist;']);
           
       end; % for n=1:nfile
       
   elseif kmen1==3; %  if making timeseries matrix
       
       
       
       % Load file specifying stns
       [file4,path4]=uigetfile('regstn*.mat','Infile specifying stns in region');
       pf4=[path4 file4];
       eval(['load ' pf4 ' Ireg;']);
       
       %    Build output filename
       fntemp = ['day' num2str(file4(7)) 'a.mat'];
       [file5,path5]=uiputfile(fntemp,'Outfile for regional tsm');
       pf5=[path5 file5];
       
       % Get outside year for max possible needed row size of daily tsm
       cgo = min(str2num(C(Ireg,50:53)));
       csp = max(str2num(C(Ireg,55:59)));
              
       % Allocate for tsm, T, and code matrix, W
       mT= (csp-cgo+1)*366; % number of rows, before trimming NaN
       nstn= length(Ireg); % number of stations
       nT= 4+ nstn; % number of cols
       T=repmat(NaN,mT,nT);
       W=T;
       nyrT = csp-cgo+1;
       
       % Initialize first 4 cols of T
       a1 = cgo:csp;
       a2=repmat(a1,366,1);
       a3=a2(:);
       T(:,1)=a3;
       b1=[(1:31)'; (1:29)'; (1:31)'; (1:30)'; (1:31)'; (1:30)'; (1:31)'; (1:31)'; (1:30)'; (1:31)'; (1:30)'; (1:31)'];
       b2=repmat(b1,nyrT,1);
       T(:,3)=b2;
       c1a=[repmat(1,31,1); repmat(2,29,1); repmat(3,31,1); repmat(4,30,1); repmat(5,31,1); repmat(6,30,1)];
       c1b=[repmat(7,31,1); repmat(8,31,1); repmat(9,30,1); repmat(10,31,1); repmat(11,30,1); repmat(12,31,1)];
       c2=[c1a; c1b];
       T(:,2)=repmat(c2,nyrT,1);
       T(:,4)=repmat((1:366)',nyrT,1);
       W(:,1:4)=T(:,1:4);
       
       
       % stninfo initialize
       stninf=[];
       
       
       for n =1:nfile; % Loop over stations
           
           % Check header against survey file info
           c=C(n,:); % char info this stn
           disp(c);
           cnm = c(9:30); % most of stn name
           cdata=fgetl(fid1); % get line from data file
           if ~strcmp(cnm,cdata(9:30));
               disp(cnm);
               disp(cdata(9:30));
               error('Miss-match on above cnm & cdata line');
           end;
           
           
           %    Allocate space for daily tsm segment
           
           mline=num1(n); % number of lines of data to read
           igo = ihead(n)+1; % first data line
           isp = igo + mline-1; % last ine
           X=repmat(NaN,31*mline,6); % more than enough space
           
           Lthis = n==Ireg;
           if any(n==Ireg); % if this is one of the stations in region
               ithis=find(Lthis);
               
               stninf=char(stninf,c);
                              
               % Loop over data lines
               
               ion=1; % initialize target line counter
               for m = 1:mline; % loop over lines for this station
                   if m==1; 
                       kopt=1;
                       yrold=[];
                       mold=[];
                   else;
                       kopt=2;
                   end;
                   
                   cline=fgetl(fid1); % get line from data file
                   [y,ycode,yrnew,mnew]=dayorg1(cline,kopt,yrold,mold);
                   
                   nvals = length(y);
                   
                   vyear = repmat(yrnew,nvals,1);
                   vmon = repmat(mnew,nvals,1);
                   vday = (1:nvals)';
                   
                   ioff=ion+nvals-1;
                   X(ion:ioff,[1 2 3  5 6])=[vyear vmon vday y ycode];
                   
                   ion=ioff+1;
                   yrold=yrnew;
                   mold=mnew;
                   
                   
                   
               end; % for m=1:mline
               
               
               %--- Trim trailing NaN rows off X
               
               L=  (all((isnan(X))'))';
               X(L,:)=[];
               [mX,nX]=size(X);
               
               
               %--- Add sequential day index (1:366)
               
               if rem(mX,366)~=0;
                   error('Row size of X not evenly divisible by 366');
               else;
                   n366=mX/366;
                   v366=repmat((1:366)',n366,1); 
                   X(:,4)=v366;
               end;
               
               
               %-- Compute target location for X precip data in T
               
               if ~[X(1,2)==1 & X(1,3)==1];
                   error(['First daily data for station ' cnm ' not Jan 1']);
               end;
               i1 = 1 +  (X(1,1)-cgo)*366; % row of T for first data
               i2 = i1 + mX - 1; % ... last data
               
               
               T(i1:i2,(ithis+4))=X(:,5); % plug in data
               W(i1:i2,(ithis+4))=X(:,6); % plug in est code
              
               
              
               
           else; % ~any(n==Ireg);  if this is not one of the stations in region
               for m = 1:mline; % loop over lines for this station
                   cline=fgetl(fid1); % get line from data file
               end;
               
               
           end; % if a station in region
           
           
           
       end; % for n=1:nfile
       
       
       % Check that no internal-row all-NaNs; if not, trim
       T1 =T(:,5:nT);
       T2=ones(mT,1);
       yr = T(:,1);
       month=T(:,2);
       day=T(:,3);
       Leap=leapyr(yr);
       
       L1 = ((all(isnan(T1')))'); % cv,1 if all values in mtx NaN this obs
       L2 = month==2 & day ==29; % cv, 1 if this Feb 29
       
       if any(L1); % if any days all-missing data
           T2(L1)=NaN;
           Ln = intnan(T2); % any of those days "internal?"
           if any(Ln);
               if Ln==1  &  ~(L2 &   ~Leap);  % internal all-missing and not Feb 29 of a non-leap year
                   error('Invalid internal NaN in T');
               else;
                   T(L1 & ~Ln,:)=[]; % remove all rows with all NaN, except "internal" rows
                   W(L1 & ~Ln,:)=[]; % remove all rows with all NaN, except "internal" rows
                   
               end;
           else;
               T(L1,:)=[];
               W(L1,:)=[];
           end;
           [mT,nT]=size(T);
       else;
       end;
         
               
       
       % Compute information for T
       I=[T(1,[1 2 3]); T(mT,[1 2 3])];
       vlist = ['T daily eischeid multi station time series matrix'];
       vlist=char(vlist,'  col 1= year');
       vlist=char(vlist,'  col 2= month');
       vlist=char(vlist,'  col 3= day of month');
       vlist=char(vlist,'  col 4= 1:366 sequential day');
       vlist=char(vlist,'  col 5= daily precip in hundredths of inches, stn 1');
       vlist=char(vlist,'  col 6= daily precip in hundredths of inches, stn 2, etc');
       vlist=char(vlist,'W daily estimation code matrix for T');
       vlist=char(vlist,'I  row 1: first year, month, day');
       vlist=char(vlist,'   row 2: last year ,month, day');
       vlist=char(vlist,'stninf stations in region');
       vlist=char(vlist,'Ireg pointer to station in station list file (big one -- Psurvey)');
              
       stninf(1,:)=[]; % trunc blank line
       % Save in .mat file
       eval(['save ' pf5 ' T W I Ireg stninf vlist;']);
       
       
       
   end; %if kmen==2; % if making individual station mat files
   
   fclose (fid1);
   
end; % switch



disp(['Number of lines = ' num2str(mX)]);
disp(['Number of files = ' num2str(nfile)]);

disp(['Elapsed time = ' num2str(etime(clock,t0)) ' sec']);

