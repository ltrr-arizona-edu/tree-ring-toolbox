% dot2crn3.m   convert Matt Therrell weird format Calif oak to crn files
%
% Matt Therrell sent data as .tab files, one per chron or merged chron, with year, sample size,
% standard chron, etc as columns.  I took that data into wordpad, copied the site info into
% mattinfo.txt, and lopped off the header lines.  The remainder of the data is all-numeric and
% can be loaded by dot2crn3.m.  Run the program separately for the std, res, and if desired, the ars
% versions of the chronologies.  Out .crn files codes with 6-letter names.  Chars `1-3 are site code.
% 4,5,6 are data type (e.g., W==width) ring part (e.g., L==latewood), and chron version (e.g., S=standard).
% Thus finwts.crn is Finley Lake, total width chron , standard

dir1='c:\data\icrns\stahle\oakresids\'; % input directory with ringwidth index spreadsheet
version='0';

kmen3=menu('Choose one',...
   'Width',...
   'Density');
switch kmen3;
case 1; % width 
   type='Width';
case 2; % density
    error('Jack, you have no density')'
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
   error('Jack, calif oaks are total rw');
case 3; % early
   width='Early';
   error('Jack, calif oaks are total rw');
otherwise;
end;


kmen2=menu('Choose one',...
   'Standard chronology',...
   'Residual chronology',...
   'ARSTAN chronology');
switch kmen2;
case 1; % standard
   ctype='Standard';
   ccode='s';
case 2; % residual
   ctype='Residual';
   ccode='r';
   
case 3; % ARSTAN
   ctype='ARSTAN';
   ccode='a';
   error('jack, you do not use arstan version');
otherwise;
end;




%---ASSEMBLE HEADER INFO

[file4,path4]=uigetfile([dir1 'infostd.txt'],'Siteinfo file');
pf4=[path4 file4];
fileinf = textread(pf4,'%s','delimiter','\n','whitespace',''); % cell matrix
Si=fileinf;
s=char(Si);
[mSi,nSi]=size(Si);

% 3 and 6-Letter code
Scode3=s(:,36:38);
s1=[type(1) width(1) version];
Scode6 = [Scode3 repmat(s1,mSi,1)];

% Site name
Sname=s(:,1:30);

% State or country
Sstate=upper(s(:,32:33));

%Lat and long
latD=s(:,46:47);
latM=s(:,49:50);
lonD=s(:,56:58);
lonM=s(:,60:61);


% Elev (m)
el = s(:,67:69);

% Collector
Guy = upper(s(:,83:97));

% Species
Species = s(:,41:44);


%***************


% Standard chrons
[file1,path1]=uigetfile([dir1 'castdcrn.dat'],'Input tsm');
pf1=[path1 file1];

file = textread(pf1,'%s','delimiter','\n','whitespace',''); % cell matrix
S=strrep(file,'    .','  NaN');
X=str2num(char(S));

% Compute number of files
yrX=X(:,1); X(:,1)=[];
[mX,nX]=size(X);
nfiles=nX; clear nX;
nsize1=mX; clear mX;

% Compute first and end year or each chron
YRon = repmat(NaN,nfiles,1);
YRoff=repmat(NaN,nfiles,1);
for n = 1:nfiles;
   x=X(:,n);
   Lgood = ~isnan(x);
   igood=find(Lgood);
   yron = yrX(min(igood));
   yroff = yrX(max(igood));
   YRon(n)=yron;
   YRoff(n)=yroff;
end;
Son = num2str(YRon,'%4.0f');
Soff=num2str(YRoff,'%4.0f');




% Loop over chronologies
for m=1:nfiles;
   
      
   
   % BUILD 3-LINE HEADER
   h1=blanks(80); h2=blanks(80); h3=blanks(80);
   h1(1:6)=Scode6(m,:);     h2(1:6)=Scode6(m,:);     h3(1:6)=Scode6(m,:);
   h1(7:8)=' 1';        h2(7:8)=' 2';        h3(7:8)=' 3';
   
   % Build output .crn file name
   sixpack=h1(1:6);
   file6=[sixpack ccode '.crn'];
   pf6=[path1 file6];
   
   species=Species(m,:);

   
   sname=upper(Sname(m,:));
   len1=length(sname);
   h1(10:(10+len1-1))=sname;
   h1(62:65)=species;
   sstate=Sstate(m,:);
   len1=length(sstate);
   h2(10:(10+len1-1))=sstate;
   elthis=el(m,:);
   h2(41:45)=[' ' elthis 'M'];
   latd=latD(m,:);
   latm=latM(m,:);
   lond=lonD(m,:);
   if strcmp(lond(1),' ');
      lond(1)='0';
   end;
   lonm=lonM(m,:);
   h2(48:49)=latd;
   h2(50:51)=[latm];
   h2(52:55)=['-' lond];
   h2(56:57)=lonm;
   yron=sprintf('%5.0f',YRon(m,:)); 
   yroff=sprintf('%5.0f',YRoff(m,:)); 
   h2(67:71)=yron;
   h2(72:76)=yroff;
   
   
   
   guy=Guy(m,:);
   len1=length(guy);
   h3(10:(10+len1-1))=guy;
   switch ctype;
   case 'Standard';
      h3(63)='_';
   case 'Residual';
      h3(63)='R';
   case 'ARSTAN';
      h3(63)='A';
   otherwise;
   end;
   
   switch type;
   case 'Width';
       switch width;
       case 'Total';
           h3(62)=' ';
       case 'Late';
           h3(62)='L';
       case 'Early'; 
           h3(62)='E';
       otherwise;
           error('in switch type');
       end;
       
   otherwise;
       error('Only width coded for so far');
   end;
      
   id=sixpack';

   % Get the series and year vector
   x=X(:,m);
   L=~isnan(x);
   yr=yrX(L);
   x=x(L);
   
   % Dummy sample size
   n=repmat(0,length(x),1);
   
   % Open file for writing
   pf2=[dir1 'stda1.crn'];
   fid=fopen(pf6,'w');
   
   
   % Header Lines
   fprintf(fid,'%s\n',h1);
   fprintf(fid,'%s\n',h2);
   fprintf(fid,'%s\n',h3);
   
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
   if rem(yr(length(yr)),10)~=0,
      rapnd=ones(9-rem(yr(length(yr)),10),1)*9990;
   end
   rfcount=rcount+length(rapnd);

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
      fprintf(fid,'%4d%3d',round(x(i)*1),n(i));
   end
   
   % Append 9990's at the end of the data if necessary
   rrem=rem(rfcount,10);
   if rrem~=0,
      for i=1:rrem,
         fprintf(fid,'%4d%3d',9990,0);
      end
   end

   
   
   
   fclose (fid);
   
end;



disp('here');


