function tempor1(hokey);
%
% Dedicated function (for moving1.m) to sort data matrix and make table of 
% worst n-year windowed periods
%
%****************  IN
% 
% D (mD x nD)r    data values, mD moving periods and nD observations
% wsize1(1 x ?)r  window widths

D=hokey{1};
YRS=hokey{2};
wsize1=hokey{3};
textin=hokey{4};
nseas=hokey{5};

namser=textin{1};
dtype=textin{2};
units=textin{3};

[mD,nD]=size(D);

% Rank D along rows
D=D';
[S,I]=sort(D);  % sorted ascending
%S=flipud(S); % sorted descending
%I=flipud(I);  % ditto


yr2=(YRS(1,1):YRS(1,2))';

% Compute year, month corresponding to index I
% Convert sequential month of end of runs to year,month
Iend=I;  % ending sequential obs of period
monend=rem(Iend,nseas);
L1temp=monend==0;
yearend=yr2(1)+floor(Iend/nseas);
monend(L1temp)=nseas;
yearend(L1temp)=yr2(1)+Iend(L1temp)/nseas-1;


ww=(wsize1(1):wsize1(2):wsize1(3));


%---- Table of Lowest value for each window size

[file2,path2]=uiputfile('*.txt','File for summary of moving window analysis');
pf2=[path2 file2];
fid2=fopen(pf2,'w');

% Loop over windows, from narrowest to widest

str1=['Series: ' namser];
str2=['Data Type: ' units];

fprintf(fid2,'%s\n\n','SUMMARY OF LOWEST VALUE FOR EACH WINDOW SIZE')
head2a='Window    End                 ';
head2b='Width     Year          Value  ';
fprintf(fid2,'%s\n',head2a);
fprintf(fid2,'%s\n\n',head2b);



fprintf(fid2,'%s\n%s\n',str1,str2);
for n =1:mD;
   howwide=ww(n);
   str11=sprintf('%2.0f',n);
   str12=sprintf('%3.0f\t     ',howwide);
   str13=sprintf('%2.0f/',monend(1,n));
   str14=sprintf('%4.0f\t',yearend(1,n));
   str15=sprintf('%g',S(1,n));
   
   if nseas>1;
      strall1=[ str12  str13 str14  str15];
   else
      strall1=[ str12  str14   str15];
   end
   
   fprintf(fid2,'%s\n',strall1);
   
end


%------------  Table of highest values

fprintf(fid2,'%s\n',blanks(10));
fprintf(fid2,'%s\n',blanks(10));
fprintf(fid2,'%s\n',blanks(10));
fprintf(fid2,'%s\n',blanks(10));
fprintf(fid2,'%s\n',blanks(10));
fprintf(fid2,'%s\n',blanks(10));

fprintf(fid2,'%s\n\n','SUMMARY OF HIGHEST VALUE FOR EACH WINDOW SIZE')
head2a='Window    End                 ';
head2b='Width     Year          Value  ';
fprintf(fid2,'%s\n',head2a);
fprintf(fid2,'%s\n\n',head2b);

for n =1:mD;
   howwide=ww(n);
   str11=sprintf('%2.0f',n);
   str12=sprintf('%3.0f\t     ',howwide);
   str13=sprintf('%2.0f/',monend(nD,n));
   str14=sprintf('%4.0f\t',yearend(nD,n));
   str15=sprintf('%g',S(nD,n));
   
   if nseas>1;
      strall1=[ str12  str13 str14  str15];
   else
      strall1=[ str12  str14   str15];
   end
   
   fprintf(fid2,'%s\n',strall1);
   
end


fclose(fid2);
