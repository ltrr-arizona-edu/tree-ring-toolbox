function rwltrim
% rwltrim:  trim .rwl file to exclude undesirable beginning or ending parts of some series
% rwltrim;
% Last revised 4-4-01
%
% You ran COFECHA and found that dating for beginning or ending of some ring-width series suspect.  You 
% do not want these questionable segments in the .rwl used to generate the chronology and sent to the 
% ITRDB.  You want to edit off the undesirable parts.
%
%*** INPUT
%
% No arguments
% You are prompted for names of input and output .rwl files
% You are prompted for name of file with list of .rw filenames, periods and comments <rwltrim?.txt>
%
%
%*** OUTPUT
%
% No arguments
% User prompted name of output .rwl fil
%
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% Easy way to build list of filenames is to edit output of DPL/SUR run on the full .rwl file
% Sometimes rwltrunc.m run after rwltrim.m to unersally lop off start years before a specified decade

[file1,path1]=uigetfile('*.rwl','Infile to be trimmed');
pf1=[path1 file1];

[file2,path2]=uiputfile('*.rwl','Outfile to hold trimmed data');
pf2=[path2 file2];

[file3,path3]=uigetfile('rwltrim?.txt','Infile with list of series & periods & comments');
pf3=[path3 file3];

D = textread(pf1,'%s','delimiter','\n','whitespace','');
D=char(D);
[mD,nD]=size(D); % mD is # rows in input .rwl
E=blanks(nD);
E=repmat(E,mD,1); % blank char matrix, same size as input .rwl
T = textread(pf3,'%s','delimiter','\n','whitespace','');

[nser,nT]=size(T); % nser is number of desired series in output .rwl

nrowout = 0; % initialize counter for number of rows in output .rwl

fid1= fopen(pf2,'w');


% Loop over series
for n = 1:nser;
    c=T{n};  % char var with sequenc no., series id, period for output .rwl, of input .rwl, and comments
    c=[blanks(1) c];
    
    lenc=length(c);
    i1=isspace(c);
    i2=diff(i1);
    i3=find(i2==-1);
    
    seqno = strtok(c(i3(1):lenc),' '); % sequence no
    cname = strtok(c(i3(2):lenc),' '); % series name
    ncname=length(cname);
    cyrgo = strtok(c(i3(3):lenc),' ');
    yrgo = str2num(cyrgo);  % start year
    cyrsp = strtok(c(i3(4):lenc),' ');
    yrsp = str2num(cyrsp); % end year
    
      
    
    % Pull block of input .rwl for this series
    L1 = strncmp(cname,cellstr(D),ncname);
    nline1=sum(L1);
    if nline1>0;
        D1=D(L1,:);
    else;
        error(['No data in ' pf1 ' for ' cname]);
    end;
    
    % Pull column of year values from rwl block. igo is the column the year begins in.  
    % Special case is like -1050, which has - in col 8.  Usually year in 9-12
    if any(D1(:,8)=='-');
        igo = 8;
    else;
        igo=9;
    end;
    isp = 12;
    yr1 = str2num(D1(:,igo:isp));
    
    
    % Check for out of range years
    if yrgo < yr1(1);
        error(['specified yrgo is ' num2str(yrgo) ' but first row of series ' cname ' is ' num2str(yr1(1))]);
    end;
    if yrsp >  (yr1(length(yr1))+9);
        error(['specified yrsp is ' num2str(yrsp) ' but last row of series ' cname ' is ' num2str(yr1(length(yr1)))]);
    end;
    
    
    c999=D1(:,16:18);  % 999 here if last row for .rwl has just a "999" (no data following)
    
    
    
    
    % Lop off any trailing rows
    L2=yr1>yrsp ;
    if any(L2);
        D1(L2,:)=[];
        yr2=str2num(D1(:,igo:isp));
    else;
        yr2=yr1;
    end;
    clear L2;
    
   
    % Lop off any leading rows 
    yrgoa = 10* floor(yrgo/10);
    L2 = yr2<yrgoa;
    if any(L2);
        D1(L2,:)=[];
        yr3=str2num(D1(:,igo:isp));
    else;
        yr3=yr2;
    end;
    clear yrgoa;
   
    
    % Fix first record if indicated
    dfirst=D1(1,:);
    if yr3(1) < yrgo;
        ndrop = yrgo-yr3(1);
        
        yrstr=num2str(yrgo,'%5.0f');
        nsize=length(yrstr);
        if igo==9;
            dfirst(igo:isp)=blanks(4);
            
        else;
            dfirst(igo:isp)=blanks(5);
        end;
        dfirst((isp-nsize+1):isp)=yrstr;
        i1 = 13; % first char to be cut
        i2 = (i1-1) + 6*ndrop; % last char to be cut
        dfirst(i1:i2)=[];
    elseif yr3(1)>yrgo;
        error(['Year on first record greater than specified yrgo for ' cname]);
    else; % Start not truncated .  no action needed
    end;
    
    if yr3(1) ~= yrgo;
        if length(dfirst)<nD;
            npad=nD-length(dfirst);
            dfirst = [dfirst blanks(npad)];
            D1(1,:)=dfirst;
        else;
            error('Should not get here');
        end;
    end;
    
    
    
    %--- Fix last record if needed
    
   % Pull last line of current version 
   [mD1,nD1]=size(D1);
   dlast=D1(mD1,:);
   
   yrspa = 10* floor(yrsp/10); % decade listed at left on line
   
   % number of trailing values to drop
   if yrgo>=0;
       ndrop = 10 - (yrsp-yrspa) -1; 
   else; % if ending year negative
       ndrop = 10 - (yrsp-yrspa) -1;
   end;
   
   if ndrop==0; % do not drop any years off last row, but must build a '999' row;
       dadd = blanks(nD);
       yradd = yrspa+10;
       stradd = num2str(yradd,'%5.0f');
       nadd = length(stradd);
       if igo==9;
           dadd(igo:isp)=blanks(4);
           dadd(1:8)=dlast(1:8);
           
       else;
           dadd(igo:isp)=blanks(5);
           dadd(1:7)=dlast(1:7);
       end;
       dadd((isp-nadd+1):isp)=stradd;
       dadd(13:18)='   999';
       D1=char(D1,dadd);
   else;
       ioff=nD; % last char to be cut
       ion = ioff - 6*ndrop +1; % first to be cut
       dlast(ion:ioff)=[];
       dlast = [dlast '   999'];
       npad = nD - length(dlast);
       if npad>0;
           dlast = [dlast blanks(npad)];
       else;
       end;
       D1(mD1,:)=dlast;
           
   end;
   [m1,n1]=size(D1);
   for m = 1:m1;
       s=deblank(D1(m,:));
       fprintf(fid1,'%s\n',s);
   end;
   
   disp(n)
    if n==12;
        disp('here');
    end;
   disp(['Finished trimming ' cname]);
end;
fclose(fid1);

