function cas2crn1(fnlist,crninf,datatype)
% cas2crn1:  casewise Stahle chronology file to .crn files
% CALL: cas2crn1(fnlist,crninf,datatype);
%
% Meko 12-15-97
%
%****************  IN *********************************************
%
% fnlist(1 x ?)s  path\filename to file with number of files in line
%     1 and paths\filenames of casewise files following. Example
%     2
%     c:\data\treedata\stahle\amrcaqdac.cas
%     c:\data\treedata\stahle\pppcadrpc.cas
%
% crninf (1 x ?)s  path\filename to file holding site information, 
%     which is assumed to be in Stahle's format
% datatype (1 x 1)s   'S', 'R', or 'A'  -- whether standard, resid, 
%     or arstan chron
%
%*******************  OUT ****************************************
%
% No output args
%
% .crn files are built, with 3-line headers.  Names of .crn files
% are automatically built from prefixes of .cas input files.  Takes
% first 7 chars and adds S, A, or R, depending on data type
% .crn files go into the current directory
%
%*******************  NOTES ****************************
%
% First worked this out for getting .crn files from Stahle's oaks.
% No sample size info is in the casewise file; zero is substituted for
% number of samples in the .crn files
%
% User runs this function once to get the standard .crn files, then again
% if he wants the other files
%
%
%---- STEPS
%
% Read initial info
% Loop over the sites
%    Compute year info
%    Get data and organize it
%    Write header lines
%    Write decade-line data
%    Save .crn fileInitial read to count number of lines in casewise file and to put NaNs
% End repeat

if strcmp(datatype,'S');
   datcode='  ';
   icol =1; % column holding data
elseif strcmp(datatype,'R');
   datcode='_R';
   icol = 2;
elseif strcmp(datatype,'A');
   datcode='_A';
   icol = 3;
else
   error('Invalid datatype arg');
end


fid7=fopen(fnlist,'r'); % open file with list of casewise filenames
fid5=fopen(crninf,'r'); % open file with long-lat, etc for files

c = fgetl(fid7);
clen=length(c);
numfiles = str2num(c);

for m = 1:numfiles;
   
    
   
   c7 = fgetl(fid7);
   pf1 = strtok(c7);
   fid1 = fopen(pf1,'r');
   
   nhead=0;
   k1=1; 
   nlines=0;
   while k1;
      c=fgetl(fid1);
      if (~feof(fid1)) | (feof(fid1) & length(c)>4);
         nlines=nlines+1;
         cstart = c(1:4);
         if strcmp(cstart,'Year');
            nhead=1;
         else
         end
      else
         k1=0;
         fclose(fid1);
         break;
      end
   end
     
  
   %-------- SIZE STORAGE 

   fid1=fopen(pf1,'r');

   X = repmat(NaN,nlines-nhead,3);
   yr = repmat(NaN,nlines-nhead,1);

   % Skip header line
   frewind(fid1);
   if nhead==1;
      c=fgetl(fid1);
      nlines=nlines-1;
   end

   for n = 1:nlines;
      c=fgetl(fid1);
      yr(n)=str2num(c(1:4));
      s = c(6:13);
      r = c(15:22);
      a = c(24:31);
      if ~isempty(str2num(s));
         X(n,1)=1000*str2num(s);
      else
         X(n,1)=NaN;
      end
      if ~isempty(str2num(r));
         X(n,2)=1000*str2num(r);
      else
         X(n,2)=NaN;
      end
      if ~isempty(str2num(a));
         X(n,3)=1000*str2num(a);
      else
         X(n,3)=NaN;
      end
   end

   %-------- Calulate start year of the three types of series

   % Check last year to see that all types have data then
   [mX,nX]=size(X);
   if any(isnan(X(mX,:)));
      error('All types of series do not have data for last year');
   end


   yrgo = yr(1);
   yrsp = yr(length(yr));
   L1 = isnan(X);
   sum1=sum(L1);
   yrgos = yrgo+sum1(1);
   yrgor = yrgo+sum1(2);
   yrgoa = yrgo+sum1(3);


   %----------- BUILD THE  .CRN FILE
   
      
   x = X(:,icol);
   x(isnan(x))=[];

   % Start and end years of data
   if strcmp(datatype,'S');
      yr1= yrgos;
   elseif strcmp(datatype,'R');
      yr1 = yrgor;
   elseif strcmp(datatype,'A');
      yr1 = yrgoa;
   else;
      error('Invalid datatype in put arg');
   end
   
   yr2=yrsp;
   nyrs = yrsp-yrgos+1;
   yr = (yr1:yr2)';

   % Compute expected number of "decade lines" of data.  Note special case
   % handling
   % of negative years.  These occur when series are BC  and the file does not
   % use the "add 8000" convention
   ii1=yr1-rem(yr1,10); % first 'decade'
   % Handle the 'negative' starting year case (a la El Malpais)
   if rem(yr1,10)<0;
      ii1=ii1-10; % start year in, say, -136, gives first decade as -140, not -130
   end
   ii2 = yr2 - rem(yr2,10);  % ending decade
   % Handle 'negative' end year
   if rem(yr2,10)<0;
      ii2=ii2-10;
   end
   nlines = round(ii2/10)-round(ii1/10) + 1;  % number of 'decade lines' of data


   % Compute cv of year labels for records
   decadd = ((0:(nlines-1))*10)';
   decbase = repmat(ii1,nlines,1);
   yrlab = decbase+decadd;
   yrlab(1)=yr1;


   % Compute number of missing values to put on first decade-line.
   % Will differ depending on whether starting year is positive or negative.
   % Negative first year example is El Malpais, NM
   if yr1>=0;
      nskip=rem(yr1,10);
      % Let's say data starts in 1652. Code would say skip first 2 values on
      %"1650" line
   else; % negative first year
      nskip = 10 - abs(rem(yr1,10)); 
      % if yrgo is -139, this says skip one value on the "-140" decade line
   end
   
   % Compute number of values to put on last decade line
   % Again, different treatment if yrsp is negative
   if yr2>=0;
      nlast = rem(yr2,10)+1;  
      ntrail = 10-nlast; % this number of 9990's
      % For example, if yr2 is 1686, says read 7 values -- 1680 to 1686
   else;  % yr2 is negative
      nlast = 10 - abs(rem(yr2,10)) + 1;
      ntrail = 10-nlast;
      % So if yr2 is -59, says read 2 values (-60 and -59) off last line
      % Or if yr2 is -51, says read 10 values (-60,-59,...,-51) of last line
   end


   % Pad series with 9990 as needed to handle start and end years
   if nskip>0;
      x= [repmat(9990,nskip,1); x];
      yr = [((yr(1)-nskip):(yr(1)-1))'; yr];
   end

   if ntrail>0;
      x = [x;  repmat(9990,ntrail,1)];
      yr = [yr; ((yr2+1):(yr2+ntrail))'];
   end





   % Write .crn file in current directory

   % Build file name
   fn1 = strtok(pf1,'.'); % cut off suffix
   fn1 = fliplr(fn1);
   fn1 = strtok(fn1,'\');
   fn1=fliplr(fn1);
   len1 = length(fn1);
   fn2=fn1;
   if len1>7;
      fn2 = fn2(1:7);
   end
   fn2 = [fn2 datatype '.crn'];

   % Build line text label
   linelab = blanks(6);
   if len1<=6;
      linelab(1:len1)=fn1(1:len1);
   else
      linelab=fn1(1:6);
   end


   % Build header lines
   c = fgetl(fid5); % read line for first file
   code6=upper(c(1:6));
   
   name=blanks(51);
   temp = strtok(upper(c(16:43)),',');
   lentemp = length(temp);
   lenuse = min(lentemp,51);
   name(1:lenuse)=temp(1:lenuse);
      
   state = blanks(13);
   temp='CALIFORNIA';
   lentemp = length(temp);
   lenuse=min(lentemp,13);
   state(1:lenuse)=temp(1:lenuse);
   
   code3=c(49:51);
   
   theman=blanks(51);
   temp='D. STAHLE';
   lentemp = length(temp);
   lenuse=min(lentemp,51);
   theman(1:lenuse)=temp(1:lenuse);

   species=blanks(16);
   temp = c(57:66);
   if strcmp(temp,'Blue Oak  ');
      specode = ' QUDO ';
   elseif strcmp(temp,'Valley Oak');
      specode = ' QULO ';
   end
   temp=upper(temp); % common name
   lentemp = length(temp);
   lenuse=min(lentemp,16);
   species(1:lenuse)=temp(1:lenuse);
      
   c = fgetl(fid5);
   latdeg = c(5:6);
   latmin = c(8:9);
   latsec = c(11:12);
   latlet = c (14); 
   if isempty(str2num(latsec));
      % no action
   else
      c1 = str2num(latdeg);
      c2 = str2num(latmin);
      c3 = str2num(latsec);
      xc = c1 +c2/60 +c3/3600;
      latdeg = sprintf('%2.0f',floor(xc));
      latmin = sprintf('%2.0f',(xc-floor(xc))*60);
   end


   londeg = c(23:25);
   lonmin = c(27:28);
   lonsec = c(30:31);
   lonlet=c(33);
   if isempty(str2num(lonsec));
      % no action
   else
      c1 = str2num(londeg);
      c2 = str2num(lonmin);
      c3 = str2num(lonsec);
      xc = c1 +c2/60 +c3/3600;
      londeg = sprintf('%3.0f',floor(xc));
      lonmin = sprintf('%2.0f',(xc-floor(xc))*60);
   end
   if strcmp(lonlet,'W');
      lonsign='-';
   else
      lonsign='+';
   end

   elevm = c(40:43);

   str1a = [code6 ' 1 '];
   str1b = sprintf('%51s',name);
   str1c = specode;
   str1d = sprintf('%s',blanks(14));
   str1 = [str1a str1b str1c str1d];

   str2a = [code6 ' 2 '];
   str2b = sprintf('%13s',state);
   str2c = sprintf('%16s',species);
   str2d = sprintf('  %4.0fM  ',str2num(elevm));
   str2e = [latdeg latmin lonsign londeg lonmin '          '];
   str2f = sprintf('%4.0f %4.0f    ',yr1,yr2);
   str2=[str2a str2b str2c str2d str2e str2f];

   str3a = [code6 ' 3 '];
   str3b = sprintf('%51s ',theman);
   str3c = datcode;
   str3d = blanks(17);
   str3=[str3a str3b str3c str3d];

   % Open file
   fid3 = fopen(fn2,'w');

   % Put the header lines on 
   fprintf(fid3,'%s\n',str1);
   fprintf(fid3,'%s\n',str2);
   fprintf(fid3,'%s\n',str3);


   % Format
   fmt1 = '%6s%4.0f';
   fmt2 = '%4.0f%3.0f';
   fmt3 = '%4.0f%3.0f\n';

   % Write Lines
   for n = 1:nlines;
      yrline = yrlab(n);
      igo = 1 + 10*(n-1);
      isp = igo+9;
      y = (x(igo:isp))';  % a rv
      ss = zeros(1,10);  % sample size (unknown);
      z = [y; ss];
      
      fprintf(fid3,fmt1,linelab,yrline);
      fprintf(fid3,fmt2,z(:,1:9));
      fprintf(fid3,fmt3,z(:,10));
   end
   fclose (fid3);
   
end
fclose all;
   
