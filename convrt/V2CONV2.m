function v2conv2
% v2conv2:  convert GHCN v2 monthly precip data and inventory files 
% CALL: v2conv2
%
% Meko 4-19-00
%
% Why?  Downloaded GHCN precip data from NCDC.  Need to create a site info 
% database and corresponding .mat files with the station data
%
%
%********************* INPUT ****************
%
% No input arguments
%
% User prompted for information on bounding lat  and lon for study area
%
%
%********************* OUTPUT ***************
%
% .MAT files with data (deg C)  in 13-col matrix X
%  .inv file with data that can be directly made into foxpro site info database
%       via excel
%
%
%******************** NOTES *****************
%
% A subset of the ghcn v2 P has adjusted records, deemed higher quality. 
% When adjusted records exist, I give these priorty in building the .mat files.
% You need to run v2conv2.m twice.  First time, generates this:
%   1. v2dbf.txt  -- the database info for adjusted stations, with the "adjusted"
%       field as 1
%   2. idadj.mat -- a string matrix of the 11-char adjusted station ids is
%      in matrix A
%   3. P1.mat, P2.mat ....    the .mat files for the adjusted stations
%
% Second time, you read the massive pv2.mat matrix and pv2.inv matrix and do the following
%   1. If a station is NOT in the adjusted set, do this:
%      - make its .mat file, with numbering sequential above the adjusted series
%      - append dbf data to v2dbf.txt, with the "adjusted" field as 0
%
% This function should be stored in the directory with the input data files. In
% matlab, user should make this the current directory.  
% 
% Have the following downloaded, unzipped, renamed files:
%   pv2.dat data file, undadjusted records
%   pv2adj.dat   data file, adjusted
%   pv2.inv, pv2adj.inv -- inventory files
%
% Output file naming convention.  I am storing the files sequentially, with
% names P1,...,P10,P11,... -- all .mat, with data in Z
%
% Output data format.  Output in the .mat files is 13 numeric columns, with
% year as col 1.  Data are in mm.
%
% I have not dealt with so-call "duplicates" yet. See readme.
%

%*********************** USER PROMPTED INFO

% Prompt for pass-1 (adjusted) or pass-2 (others)
kpass=menu('Choose One',...
   'Pass-1: adjusted stations',...
   'Pass-2: other stations');

% Trial run or the real thing?
kdebug = menu('Debugging mode or the real thing?',...
   'Debugging', 'The Real Thing');
if kdebug==1 & kpass==2;
   error('debug mode not valid for a pass-2 run');
end;

% If a pass-2 run, read the id-file from pass-1 and compute starting index for
% sequential station numbering
if kpass~=1; % this is a pass-2 run
   load idadj; % ids of adjusted station in A
   nadjust = size(A,1); % number of adjusted stations);
   istart=nadjust+1; %  starting file number for .mat file
else;
   istart=1;
end;


dtype='pcp'; % data is precip
dcode='P';  % .mat output files begin with this letter
if kpass==1; % adjusted stations only
   fndat='pv2adj.dat'; % input data file
else;
   fndat='pv2.dat'; % input data file
end;
%fndat='dopey.dat'; % input data file
pathmat='c:\data\ghcn\pcpv2\matfls\';


%---------- Long-lat boundaries of study area

prompt={'Southernmost latitude','Northernmost latitude',...
      'Easternmost longitude','Westernmost Longitude'};
% default domain
if kdebug~=1;
   deftemp={'20','70','-96','-155'}; 
else;
   if kpass==2; 
      error('Pass-2 mode invalid if debug mode');
   end;
   %deftemp={'31.0','34.5','35','-7'}; % debug-mode lon-lat domain
   deftemp={'52.2','54.3','70','61'}; % debug-mode lon-lat domain
end;

tittemp='Specify Boundary of Study Area';
lineNo=1;
answer=inputdlg(prompt,tittemp,lineNo,deftemp);
lats=str2num(answer{1}); 
latn=str2num(answer{2}); 
lone=str2num(answer{3}); 
lonw=str2num(answer{4}); 

% Cleanup
clear answer prompt tittemp deftemp lineNo kraw ktemp dtype

%**********  BUILD LONG-LAT-SCREENED INVENTORY FILE AND ASSOCIATED STRING MATRIX
%  OF 11-CHARACTER ALL-NUMERIC STATION CODES 

clc
disp('Building the screened inventory file and string matrix of station codes...');

%--------------  Open full input inventory file
if kpass==1;
   file1 = textread('pv2adj.inv','%s','delimiter','\n','whitespace','');
else;
   file1 = textread('pv2.inv','%s','delimiter','\n','whitespace','');
end;
file1=char(file1); % cell to char array

% Build logical pointer to rows of stations in lon-lat range
clat=file1(:,44:49);  clon=file1(:,51:57); % pull cols of lat and lon
cid = file1(:,1:11); % station id
celev=file1(:,59:62);
cname=file1(:,13:43); % station name
lat=str2num(clat); lon=str2num(clon); % numeric
elev=str2num(celev);
L1=lat>=lats & lat<=latn;
L2=lon>=lonw & lon<=lone;
L3=L1 & L2;
if sum(L3)<1;
   error('No stations in specified lat-lon range');
   nstn1=[];
else;
   nstn1=sum(L3);
   disp([int2str(nstn1) ' stations in lat-lon range']);
end;

% Form matrices of data for stations lon-lat range
cid=cid(L3,:);
cname=cname(L3,:);
file1=file1(L3,:);
D=[lat(L3) lon(L3) elev(L3)]; 

%  IF PASS-2, CULL STATIONS NOT IN PASS-1 SET
if kpass==2;
   Ladj=logical(zeros(nstn1,1)); % initialize logical pointer to stations also in pass-1 set
   for n=1:nstn1;
      cthis =cid(n,:);
      ilook = strmatch(cthis,A);
      if length(ilook)==1;
         Ladj(n)=1;
      elseif length(ilook)>1;
         error([cthis ' has more than one match in adjusted-set ids A']);
      else;
         % no action needed, this station not in adjusted set
      end;
   end;
   if any(Ladj);
      cid(Ladj,:)=[];
      cname(Ladj,:)=[];
      file1(Ladj,:)=[];
      D(Ladj,:)=[];
   else;
      % no action needed.  None of the full-set stations are also in pass-1 set
   end;
   if all(Ladj);
      error('Wait, all stations flagged as being in pass-1 set');
   end;
   nstn1= size(file1,1); % revised number of acceptable stations. Revised so that
   % includes stations in lat-lon range and NOT in adjusted (pass-1) set
end;

clear L1 L2 L3 clat clon celev lat lon elev Ladj cthis ilook;

%%fid1 = fopen('v2inv.txt','r'); % primary input inventory
%%fid2 = fopen('v2inv.sb1','w'); % intermediate inventory; records in lat-long range

if kpass==1;
   disp(['   Done -- ' int2str(nstn1) ' adjusted records in the lon-lat screened inventory files']);
else;
   disp(['   Done -- ' int2str(nstn1) ' additional non-adjusted records in inventory']);
end;



%*************  CULL DATA RECORDS (IF THEY EXIST) FOR SCREENED INVENTORY
%
% Approach is to read data file record by record, each time comparing the
% station id (first 11 chars) with the previously culled set of acceptabe
% station ids. A hit means
% the data record being read belongs to a desired station. So far
% I do not handle ghcn "duplicates" (see readme), but I do save the 
% character indicating whether a staton has duplicates, so that my info
% database can be searched for stations that might have duplicates.
% The 12th character in the data records indicates presence of dupes.
%
% In process, will build a third inventory file for the station info of the 
% selected stations, and a string matrix of same row size with the station codes

% Recall, that lat-lon screened inventory file is file1, with ids in cid

% Open outfile for final inventory: not only is station in long-lat search region,
% but monthly data exists
fid3 = fopen('v2inv.sb2','w'); 

fid5=fopen(fndat,'r'); % open input file of monthly data

% Open intermediate output file to hold monthly alphanumeric records for
% selected stations
fid6=fopen('v2dat.sb1','w'); 


disp('Reading data records and building first intermediate data file');
nstn2=0; % initialize counter of selected stations
S2 = repmat(blanks(12),nstn1,1); % initialize string matrix of station codes
Sline2=repmat(blanks(62),nstn1,1);

nlines=zeros(nstn1,1); % initialize to hold number of data lines in final
% selected station segments

lncount=0;
cid2old = blanks(12); % initialize 12-char id of 'previous' record

kwhile1 = 1;
while kwhile1;
   c = fgetl(fid5);
   disp(c);
   lncount=lncount+1;
   %if kdebug==1 & (lncount==100 | lncount==200 | lncount==300);
      %disp(c);
      %disp('Press any Key to continue');
      %ause
   %end
   
      
   if ~(feof(fid5) | length(c)<50 | (kdebug==1 & lncount>1000));
      cid1 = c(1:11);
      cdupe =c(12); % duplicate flag. ==1 means none, ==0 means yes
      cid2=[cid1 cdupe];
      i1 = strmatch(cid1,cid);
      nmatch=length(i1);
      if nmatch>0;
         if nmatch>1;
            fclose all;
            error('More than one id match');
         else; % Just 1 id match -- the way it should be
            
            % Write line of data
            fprintf(fid6,'%s\n',c);
            if strcmp(cid2,cid2old); % same id as previous record
               %  No need to increment station counter or store id in inventory
               nlines(nstn2)=nlines(nstn2)+1;
            else; % New id -- new station found
               nstn2 = nstn2+1;
               disp(['Found-station count = ' int2str(nstn2)]);
               cid2old=cid2; % re-set 'previous' id
               S2(nstn2,:)=cid2; % store 12-digit id
               nlines(nstn2)=1;
               
               % Get the full line of the inventory info, and write to file; store too
               Sline2(nstn2,:)=file1(i1,:);
               sline=file1(i1,:);
               fprintf(fid3,'%s\n',sline); % Write station info line to inventory file
            end
         end
      end
      
         
   else
      kwhile1=0;
   end
end

% cleanup
clear sline c cid1 cid2 cid2old cdupe nid L1 lncount ID1 ntemp kwhile1

disp(['   Done, and ' int2str(nstn2) ' stations in final inventory']);;

nlines=nlines(1:nstn2);
S2=S2(1:nstn2,:);
Sline2=Sline2(1:nstn2,:);

fclose all;


%**********  LOOP OVER THE nstn2 SELECTED STATIONS, DOING THE FOLLOWING
%
% Build filename of .mat file to hold the data
% Convert string data to numeric, and store numeric matrix for the station
% Replace missing values (-9999) with NaNs
% Replace trace (-8888)with 1's (which is 1/10 mm), since input data are 10ths of mm
% Insert any necessary all-NaN years to make rows continuous
% Convert from 10ths of mm to mm
% Compute and store the first and last years of 
%   All data
%   The period with no missing monthly data
% Save .mat file of monthly data, in X
% Build ascii database record
% Add record to output ascii database site info file
% ...............................................................


fid7 =fopen('v2dat.sb1','r');

if kpass==1;
   fid8 = fopen('v2dbf.txt','w'); 
else;
   fid8=  fopen('v2dbf.txt','a');
end;


disp(['Building the .mat files and final dbf import file']);
for n = 1:nstn2;  % Loop index will be part of .mat file name
   disp(['   Working on station ' int2str(n) ' of ' int2str(nstn2)]);
   s2=S2(n,:);  % 12-digit site + dupe code
	s3=Sline2(n,:);  % full line of inventory info
   fnmat = [dcode int2str(n+istart-1) ]; % name of output.mat file
   n1 = nlines(n);  % number of lines of data to read
   X1 = repmat(NaN,n1,12);
   yr1 = repmat(NaN,n1,1);
   for m = 1:n1; % loop over data records
      xmon = repmat(NaN,1,12); % allocate
      c = fgetl(fid7); % get a line of data
      % Check that id matches
      if ~strcmp(s2(1:12),c(1:12));
         disp(['Site ' int2str(n) ': ' s2]);
         fclose all
         error('Id mismatch for above site');
      end
      
      yr1(m) = str2num(c(13:16));
      
      % Store one month of data
      for k = 1:12; % loop over months
         igo = 17 + (k-1)*5;
         isp = igo + 4;
         cdat = c(igo:isp);
         xmon(k)=str2num(cdat);
      end
      X1(m,:)=xmon;
   end
   
   % Convert -9999 to NaN
   L1 = X1==-9999;
   if any(any(L1));
      X1(L1)=NaN;
   end
   
   % Convert -8888 to 1 (trace to 1/10 mm)
   L2 = X1==-8888;
   if any(any(L2));
      X1(L2)=1;
   end

      % Convert from tenth of mmm to mm
   X1 = 0.1 * X1;
   
   
   %--------  Insert any missing years
   yrgo1 = yr1(1); yrsp1=yr1(n1);  % first and last year with monthly data
   diff1 = diff(yr1);
   if ~all(diff1==1); % yr1 does not increment in steps of 1 year throughout
      yr2 = (yrgo1:yrsp1)';
      nyrs2 = length(yr2);
      X2 = repmat(NaN,nyrs2,12);
      iyr = yr1 - yrgo1 + 1; % index of existing data into X1
      X2(iyr,:) = X1;
   else
      X2=X1;
      yr2=yr1;
   end
   
   %---------  Find first, last year of period with no missing data
   [yrgo2,yrsp2]=nonan1(X2,yr2);   
   
   %--------------- Save .mat file of ascii data
      X=[yr2 X2];
   eval(['save ' pathmat fnmat ' X;']);
      
   
   %-------	BUIILD DATABASE FIELDS

	str1 = sprintf('%5s\t',fnmat); % .mat file name
   str2 = sprintf('%3s\t',s3(1:3)); % country code
   str3 = sprintf('%5s\t',s3(4:8)); % wmo station id
   str4 = sprintf('%3s\t',s3(9:11)); % near-to-wmo indicator
   if kpass==1;
      str4a=sprintf('%1s\t','1'); % an 'adjusted' station
   else;
      str4a=sprintf('%1s\t','0'); % unadjusted station
   end;
      str5 = sprintf('%1s\t',s2(12)); % duplicate indicator
   str6 = sprintf('%s\t',s3(13:43)); % station name
   
   str7 = sprintf('%6.2f\t',str2num(s3(44:49))); % decimal deg lat
   str8 = sprintf('%7.2f\t',str2num(s3(51:57))); % decimal deg lon
   
   % Given elevation (m), -999 if not available
   eltemp = str2num(s3(59:62));
   if eltemp==-999;
      eltemp=-999;
   end
   str9 = sprintf('%4.0f\t',eltemp);    
            
   % First and last years with any data
   stryrgo1 = sprintf('%4.0f\t',yrgo1);
   stryrsp1 = sprintf('%4.0f\t',yrsp1);
   
   % First and last years of longest continuous period with
   % no monthly data missing (most recent period if a tie)
   stryrgo2 = sprintf('%4.0f\t',yrgo2);
   stryrsp2 = sprintf('%4.0f\t',yrsp2);
          
   % Concatenate strings
   strsub1 = [str1 str2 str3 str4  str4a str5 str6 str7 str8 str9 ];
   strsub2 = [stryrgo1 stryrsp1 stryrgo2 stryrsp2];
   strall = [strsub1 strsub2 ];
   
   fprintf(fid8,'%s\n',strall);
end

   

fclose all;

