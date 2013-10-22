function v2conv1
% v2conv1:  convert GHCN v2 monthly temperature data and inventory files 
% CALL: v2conv1
%
% Meko 11-13-97
%
% Why?  Downloaded GHCN tmp data from NCDC.  Need to create a site info 
% database and corresponding .mat files with the station data
%
%
%********************* INPUT ****************
%
% No input arguments
%
% User prompted for information on bounding lat  and lon for study area,
% and for choices or 'raw' vs 'adjusted' and whether max or min data
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
% This function should be stored in the directory with the data files. In
% matlab, user should make this the current directory
% 
% Have only one downloaded inventory file, but 6 different data files, with
% various combinations of (1) raw and adjusted, and (2) max, min, and mean
% I use the following coding of .mat files:
%
% A,B,C = max, min, mean of raw
% D,E,F = max, min, mean of adjusted
%
% A00010.mat -- max tmp, station sequence 00001 in my selected set, and 
% NCDC 'duplicate' number 0 (last character).  The file A00010.mat contains
% the 13-col temperature data in matrix X.
%
%
%
% I handle only the max and min series so far.  Mean series not year
% programmed.
% With the 'mean' series, must deal with 'duplicates'.  I  **plan to** use simple approach
% of treating stations with duplicates.  I take series with most monthly data
% as starting point, then use the duplicate series to to replace any missing
% values or to extend the series. This merging will only be necessary for the
% mean series.  Apparently, the max and min series have only one time series
% per station.
%
%************************ END OF COMMENTS **************************





%*********************** USER PROMPTED INFO

% Trial run or the real thing?
kdebug = menu('Debugging mode or the real thing?',...
   'Debugging', 'The Real Thing');


% Whether raw or adjusted data
ktemp = menu('Type of temperature data','Raw','Adjusted');
if ktemp==1;
	kraw = 1;
else;
	kraw = 0;
end


% Max, min or mean
ktemp=menu('Maximum, Minimum, or Mean?','Max','Min', 'Mean');
if ktemp==1;
	dtype='max';
elseif ktemp==2;
	dtype='min';
else;
	error('v2conv1.m only handles max or min data so far -- not mean data');
	dtype='mean';
end
% Set 1-letter data type code
if kraw==1;
	if strcmp(dtype,'max');
      dcode='A';
      fndat = 'A.dat';
	elseif strcmp(dtype,'min');
      dcode='B';
      fndat = 'B.dat';
	else;
      dcode='C';
      fndat='C.dat';
	end
else;
	if strcmp(dtype,'max');
      dcode='D';
      fndat='D.dat';
	elseif strcmp(dtype,'min');
      dcode='E';
      fndat='E.dat';
	else
      dcode='F';
      fndat='F.dat';
	end
end

kmen1=menu('Type of data',...
   'Max adj tmp (mat files to winbook)',...
   'Min adj tmp (mat files to tree)');
if kmen1==1,
   pathmat = ['c:\data\ghcn\tmpv2_0\matfls' dcode '\'];
elseif kmen1==2;
   pathmat = ['o:\clim\ghcn\turbe\matfls' dcode '\'];
else;
   error('so far only use max and min adj');
end



%---------- Long-lat boundaries of study area

prompt={'Southernmost latitude','Northernmost latitude',...
	'Easternmost longitude','Westernmost Longitude'};
deftemp={'20','70','-96','-155'};
tittemp='Specify Boundary of Study Area';
lineNo=1;
answer=inputdlg(prompt,tittemp,lineNo,deftemp);
lats=str2num(answer{1}); 
latn=str2num(answer{2}); 
lone=str2num(answer{3}); 
lonw=str2num(answer{4}); 

% Override boundaries for debugging mode
if kdebug==1;
   lats=31.0; latn=34.5; lonw=-7; lone=35;
end



% Cleanup
clear answer prompt tittemp deftemp lineNo kraw ktemp dtype


%**********  BUILD LONG-LAT-SCREENED INVENTORY FILE AND ASSOCIATED STRING MATRIX
%  OF 11-CHARACTER ALL-NUMERIC STATION CODES 

clc
disp('Building the screened inventory file and string matrix of station codes...');

%--------------  Open inventory file and data file
fid1 = fopen('v2inv.txt','r'); % primary input inventory
fid2 = fopen('v2inv.sb1','w'); % intermediate inventory; records in lat-long range

blnk11=blanks(11);
blnk12=blanks(12);
S1 = repmat(blnk11,8000,1); % to store intermediate inventory station codes
Sline = repmat(blanks(100),8000,1); % to store full line of inventory info

kwhile1=1;
nstn1 = 0;
while kwhile1==1;
   c=fgetl(fid1);
   if ~(feof(fid1) | length(c)<100);
      clat = str2num(c(44:49));
      clon = str2num(c(51:57));
      if clat >= lats & clat<=latn & clon>=lonw & clon<=lone; 
         nstn1 = nstn1+1;
         fprintf(fid2,'%s\n',c);
         S1(nstn1,:)=c(1:11);
         Sline(nstn1,:)=c(1:100);
                  
      end
   else
      kwhile1=0;
   end
   
end

S1 = S1(1:nstn1,:);
Sline=Sline(1:nstn1,:);

disp(['   Done -- ' int2str(nstn1) ' records in the screened files']);
% Cleanup
clear clat clon kwhile1  c
fclose (fid1); % close initial inventory file


%*************  CULL DATA RECORDS (IF THEY EXIST) FOR SCREENED INVENTORY
%
% Approach is to read data file record by record, each time comparing the
% station id (first 11 chars) with the full set of ids in S1. A hit means
% the data record being read belongs to a desired station. So far
% I treat only max and min data. If I ever extend to mean data, will need to
% consider that because of 'retained duplicates', multiple time series could
% be there for one station.  The 12th character in the data records would be
% used to distinguish duplicates.  
% 
%
% In process, will build a third inventory file for the station info of the 
% selected stations, and a string matrix of same row size with the station codes

% Close the intermediate inventory file, then reopen for writing
fclose(fid2); 
fid2=fopen ('v2inv.sb1','r');

% Open outfile for final inventory: not only is station in long-lat search region,
% but monthly data exists
fid3 = fopen('v2inv.sb2','w'); 

fid5=fopen(fndat,'r'); % open input file of monthly data

% Open intermediate output file to hold monthly alphanumeric records for
% selected stations
fid6=fopen('v2dat.sb1','w'); 


% Make a numberic version of S1
disp('Building a numeric version of first intermediate inventory file');
ID1=zeros(nstn1,1);
for ntemp = 1:nstn1;
   c = S1(ntemp,:);
   ID1(ntemp)=str2num(c);
end
disp('   Done');


disp('Reading data records and building first intermediate data file');
nstn2=0; % initialize counter of selected stations
S2 = repmat(blnk12,nstn1,1); % initialize string matrix of station codes
Sline2=repmat(blanks(100),nstn1,1);

nlines=zeros(nstn1,1); % initialize to hold number of data lines in final
% selected station segments

lncount=0;
cid2old = blnk12; % initialize 12-char id of 'previous' record

kwhile1 = 1;
while kwhile1;
   c = fgetl(fid5);
   lncount=lncount+1;
   %if kdebug==1 & (lncount==100 | lncount==200 | lncount==300);
      %disp(c);
      %disp('Press any Key to continue');
      %ause
   %end
   
      
   if ~(feof(fid5) | length(c)<50 | (kdebug==1 & lncount>1000));
      cid1 = c(1:11);
      cdupe =c(12);
      cid2=[cid1 cdupe];
      nid = str2num(cid1);
      L1= nid==ID1;
      if any(L1);
         if sum(L1)~=1;
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
               cid2old=cid2; % re-set 'previous' id
               S2(nstn2,:)=cid2; % store 12-digit id
               nlines(nstn2)=1;
               
               % Get the full line of the inventory info, and write to file; store too
               Sline2(nstn2,:)=Sline(L1,:);
               sline=Sline(L1,:);
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



%**********  LOOP OVER THE nstn2 SELECTED STATIONS, DOING THE FOLLOWING
%
% Build filename of .mat file to hold the data
% 
% Convert string data to numeric, and store numeric matrix for the station
% Replace missing values with NaNs
% Insert any necessary all-NaN years to make rows continuous
% Compute and store the first and last years of 
%   All data
%   The period with no missing monthly data
% Save .mat file of monthly data, in X
% Build ascii database record
% Add record to output ascii database site info file
% ...............................................................


fid7 =fopen('v2dat.sb1','r');
fid8 = fopen('v2dbfd.txt','w'); 

a=NaN;

for n = 1:nstn2;
   s2=S2(n,:);  % 12-digit site + dupe code
	s3=Sline2(n,:);  % full line of inventory info
   fnmat = [dcode int2str(n) ]; % name of output.mat file
   n1 = nlines(n);  % number of lines of data to read
   X1 = repmat(a,n1,12);
   yr1 = repmat(a,n1,1);
   for m = 1:n1; % loop over data records
      xmon = repmat(a,1,12); % allocate
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
      X1(L1)=a;
   end
   
   % Convert from tenth of deg C to deg C
   X1 = 0.1 * X1;
   
   
   %--------  Insert any missing years
   yrgo1 = yr1(1); yrsp1=yr1(n1);  % first and last year with monthly data
   diff1 = diff(yr1);
   if ~all(diff1==1); % yr1 does not increment in steps of 1 year throughout
      yr2 = (yrgo1:yrsp1)';
      nyrs2 = length(yr2);
      X2 = repmat(a,nyrs2,12);
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
   str5 = sprintf('%1s\t',s2(12)); % duplicate indicator
   str6 = sprintf('%s\t',s3(13:42)); % station name
   
   str7 = sprintf('%6.2f\t',str2num(s3(44:49))); % decimal deg lat
   str8 = sprintf('%7.2f\t',str2num(s3(51:57))); % decimal deg lon
   
   % Given elevation (m), -999 if not available
   eltemp = str2num(s3(59:62));
   if eltemp==-999;
      eltemp=-999;
   end
   str9 = sprintf('%4.0f\t',eltemp);    
   
   str10 = sprintf('%4.0f\t',str2num(s3(64:67))); % interpoated el (m)
         
   % First and last years with any data
   stryrgo1 = sprintf('%4.0f\t',yrgo1);
   stryrsp1 = sprintf('%4.0f\t',yrsp1);
   
   % First and last years of longest continuous period with
   % no monthly data missing (most recent period if a tie)
   stryrgo2 = sprintf('%4.0f\t',yrgo2);
   stryrsp2 = sprintf('%4.0f\t',yrsp2);
   
   str11 = sprintf('%s\t',s3(68)); % urban code
   
   % Population (thousands), or -9 for not applicable
   pop = str2num(s3(70:73));
   if pop==-9;
      pop=-9;
   end
   str12= sprintf('%4.0f\t',pop);
   
   str13 = sprintf('%2s\t',s3(74:75)); % topo code
   str14 = sprintf('%2s\t',s3(76:77)); % veg code
   str15 = sprintf('%2s\t',s3(78:79)); % location code
   
   % Km from coast, if coastal station; otherwise -9
   dcoast = str2num(s3(80:81));
   if dcoast==-9;
      dcoast=-9;
   end
   str16=sprintf('%2.0f\t',dcoast'); 
   
   str17 = sprintf('%s\t',s3(82)); % airport(A) or not (x)
   
   % Distance airport to urban area;  -9 if not airport
   dtown = str2num(s3(83:84));
   if dtown==-9;
      dtown=-9;
   end
   str18 = sprintf('%3.0f\t',dtown);
   
   
   % Vegetation description
   str19 = sprintf('%s ', s3(85:100));
   
   % Concatenate strings
   strsub1 = [str1 str2 str3 str4 str5 str6 str7 str8 str9 str10];
   strsub2 = [stryrgo1 stryrsp1 stryrgo2 stryrsp2];
   strsub3 = [str11 str12 str13 str14 str15 str16 str17 str18 str19];
   strall = [strsub1 strsub2 strsub3];
   
   fprintf(fid8,'%s\n',strall);
end

   

fclose all;

