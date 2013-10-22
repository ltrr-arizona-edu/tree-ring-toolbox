% dot2crn2.m   convert Matt Therrell weird format latewood files to crn files
%
% Matt Therrell sent data as .tab files, one per chron or merged chron, with year, sample size,
% standard chron, etc as columns.  I took that data into wordpad, copied the site info into
% mattinfo.txt, and lopped off the header lines.  The remainder of the data is all-numeric and
% can be loaded by dot2crn2.m.  Run the program separately for the std, res, and if desired, the ars
% versions of the chronologies.  Out .crn files codes with 6-letter names.  Chars `1-3 are site code.
% 4,5,6 are data type (e.g., W==width) ring part (e.g., L==latewood), and chron version (e.g., S=standard).
% Thus MCWWLS.crn is McCurtain county  standard chron of latewood width

dir1='c:\projs\ai3\prwcleav\'; % input directory with latewood files as culled from pmail
% also used for output

kmen3=menu('Choose one',...
   'Width',...
   'Density');
switch kmen3;
case 1; % width 
   type='Width';
case 2; % density
   type='Density';
otherwise;
end;


kmen1=menu('Choose one',...
   'Total ring width',...
   'Latewood width',...
   'Earlywood width');
switch kmen1;
case 1; % total ringwidth
   width='Total';
case 2; % latewood sidth
   width='Late';
case 3; % early
   width='Early';
otherwise;
end;


kmen2=menu('Choose one',...
   'Standard chronology',...
   'Residual chronology',...
   'ARSTAN chronology');
switch kmen2;
case 1; % standard
   ctype='Standard';
case 2; % residual
   ctype='Residual';
case 3; % ARSTAN
   ctype='ARSTAN';
otherwise;
end;




%---ASSEMBLE HEADER INFO

[file4,path4]=uigetfile('mattinf.txt','Siteinfo file');
pf4=[path4 file4];
fileinf = textread(pf4,'%s','delimiter','\n','whitespace',''); % cell matrix
Si=fileinf;
s=char(Si);
[mSi,nSi]=size(Si);

% 3 and 6-Letter code
Scode3=s(:,39:41);
s1=[type(1) width(1) ' '];
Scode6 = [Scode3 repmat(s1,mSi,1)];

% Site name
Sname=s(:,1:29);

% State or country
Sstate=upper(s(:,32:37));

%Lat and long
latD=s(:,44:45);
latM=s(:,47:48);
lonD=s(:,51:53);
lonM=s(:,55:56);

% Elev (m)
el = s(:,60:63);

% Collector
Guy = upper(s(:,72:99));

% Species
Species = s(:,66:69);

% INput filenames
Infiles = [Scode3 repmat('LW.dat',mSi,1)];
nfiles=mSi;

% LOOP OVER CHRONOLOGY FILES
for n = 1:nfiles;
   fn=[dir1 Infiles(n,:)]; % file name for input
   eval(['load ' fn]); % load input data
   
   % Store data 
   f = strtok(fn,'.'); % remove suffix
   f = fliplr(f);
   f = strtok(f,'\');
   f = fliplr(f);
   eval(['Z = ' f ';']);
   
   % Note 8 cols in Z
   % 1 year
   % 2 n cores
   % 3,4 ?
   % 5 raw index
   % 6 std index
   % 7 res index
   % 8 arstan index
   
   
   % STANDARD CHRON
   
   fnout = [dir1 f ctype(1) '.crn']; % output filename
      
   % Get the series and year vector
   switch ctype;
   case 'Standard';
      x=Z(:,6)*1000;
   case 'Residual';
      x=Z(:,7)*1000;
   case 'ARSTAN';
      x=Z(:,8)*1000;
   otherwise;
   end;
   
   yr=Z(:,1);
   nc=Z(:,2);
   L=~isnan(x);
   yr=yr(L);
   x=x(L);
   nc=nc(L);
   
   % Compute first and end year 
   yron=min(yr);
   yroff=max(yr);
   son=sprintf('%5.0f',yron);
   soff=sprintf('%5.0f',yroff);
      
   % Convert all the NaN's in X into 9990's
   lmx=isnan(x);
   if any(lmx);
      x(lmx)=9.990;
   end;
   
   bcount=1;
   while lmx(bcount)==1,
      bcount=bcount+1;
   end
   bcount=bcount-1;
   rcount=1;
   k=length(lmx);
   while lmx(k)==1,
      k=k-1;
      rcount=rcount+1;
   end
   rcount=rcount-1;

   % Append 9990's in the beginning of the data if necessary
   if rem(yr(1),10)~=0,
      fapnd=ones(rem(yr(1),10),1)*9990;
   end
   bfcount=bcount+length(fapnd);
   
   % Append 9990's at the end of the data if necessary
   if rem(yr(length(yr))+1,10)~=0,
      rapnd=ones(10-rem(yr(length(yr))+1,10),1)*9990;
      rfcount=rcount+length(rapnd);
   else;
      rfcount=rcount+0;
   end
   
   
   % Open file for writing
   fid=fopen(fnout,'w');
   
   % BUILD 3-LINE HEADER
   Scode6(n,6)=upper(ctype(1));
   h1=blanks(80); h2=blanks(80); h3=blanks(80);
   h1(1:6)=Scode6(n,:);     h2(1:6)=Scode6(n,:);     h3(1:6)=Scode6(n,:);
   h1(7:8)=' 1';        h2(7:8)=' 2';        h3(7:8)=' 3';
   id=Scode6(n,:);
   
   sname=Sname(n,:);
   species=Species(n,:);
   len1=length(sname);
   h1(10:(10+len1-1))=sname;
   h1(62:65)=species;
   sstate=Sstate(n,:);
   len1=length(sstate);
   h2(10:(10+len1-1))=sstate;
   elthis=el(n,:);
   h2(41:45)=[elthis 'M'];
   latd=latD(n,:);
   latm=latM(n,:);
   lond=lonD(n,:);
   if strcmp(lond(1),' ');
      lond(1)='0';
   end;
   
   lonm=lonM(n,:);
   h2(48:49)=latd;
   h2(50:51)=[latm];
   h2(52:55)=['-' lond];
   h2(56:57)=lonm;
   h2(67:71)=son;
   h2(72:76)=soff;
   
   guy=Guy(n,:);
   len1=length(guy);
   h3(10:(10+len1-1))=guy;
   switch ctype;
   case 'Standard';
      h3(63)=' ';
   case 'Residual';
      h3(63)='R';
   case 'ARSTAN';
      h3(63)='A';
   otherwise;
   end;
      
   % Header Lines
   fprintf(fid,'%s\n',h1);
   fprintf(fid,'%s\n',h2);
   fprintf(fid,'%s\n',h3);


      
      
   % Append 9990's in the beginning of the data if necessary
   brem=rem(bfcount,10);
   if brem~=0,
      bwhln=fix(bfcount/10);
      if bwhln==0,
         fprintf(fid,'%s',id);
         fprintf(fid,'%4d',yr(1));
         for i=1:length(fapnd),
            fprintf(fid,'%4d%3d',fapnd(i),0);
         end
      else
         fprintf(fid,'%s',id);
         fprintf(fid,'%4d',yr(bwhln*10+brem-length(fapnd)+1));
         for i=1:brem,
            fprintf(fid,'%4d%3d',9990,0);
         end
      end
   end
   
   % Main loop : write the actual data x
   for i=bcount+1:length(x)-rcount,
      if rem(yr(i),10)==0,
         fprintf(fid,'\n%s%4d',id,yr(i));
      end
      fprintf(fid,'%4d%3d',round(x(i)*1),nc(i));
   end
   
   % Append 9990's at the end of the data if necessary
   rrem=rem(rfcount,10);
   if rrem~=0,
      for i=1:rrem,
         fprintf(fid,'%4d%3d',9990,0);
      end
   end

   
   fclose(fid);

   
end;
