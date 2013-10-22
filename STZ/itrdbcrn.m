function itrdbcrn
% itrdbcrn: writes tree-ring indices to an itrdb-formated .crn file
% itrdbcrn;
% Last revised 1-9-01, to handle standard and residual indices in one run
%
% Writes standard and residual tree-ring index to ITRDB-formatted .crn files. User 
% prompted for measurement type (width total, earlywood, latewood), and for whether 
% sample size should be number of cores or number of trees in any year. User  
% optionally can add ITRDB-formatted 3-line header.
%
%*** INPUT
%
% You are prompted for several pieces of information:
%   -point to the input .mat storage file with the indices
%   -a 6-character site id to use as the file name prefix
%   -whether the .crn file is to contain standard or residual index, and 
%   -whether number of cores or number of trees should be used as the sample size
%   -whether to build a 3-line "official" header -- if so, you will be
%    prompted for species, lat/lon, etc.
%
%*** OUTPUT 
% 
% Produced is an output .crn file in the format required by the ITRDB. You will
% still need to run "header" to put the header information on the file
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED -- none
%
%*** TOOLBOXES NEEDED -- none
%
%*** NOTES
%
% itrdbcrn.m is usually the last Matlab step chronology development
% The file has the following format:
%  Leftmost 6 chars are site ID string (e.g., PADSTD)
%  Next 4 columns contain "decade year"
%  Remaining columns are the data and sample size in individual years
%
% Output filename.  If the output file is the standard index, the
% file name is the xxxxxx.crn, where xxxxxx is the 6 char id prompted for. If 
% data is residual chron, filename has "r" appended (e.g., fndxxxr.crn)
%
% Output filename for earlywood and latewood not coded by any appended letter, but 
% usually within the 6-letter code.  Say xxx is the site code. The total width, 
% early width and late width files, standard version, are then xxxwt1.crn,
% xxxwe1.crn and xxxwl1.crn.
%
% You can add the ITRDB header lines (3 lines) by supplying the complete information 
% on species, lat/long etc., or by copying th 3 header lines from another .crn file 
% and optionally modifying  the site name, species, and data code.
%
% Note that col 63, line 3 of header has blank for standard chron and R for residual.
   

%-- Prompt for input .mat file
[file1,path1]=uigetfile('*.mat','Input .mat file with site index');
pf1=[path1 file1];

%-- Prompt for measurement variable
kmeas=menu('Which type of measurement?',...
   'Total Ring Width',...
   'Latewood Width',...
   'Earlywood Width');
if kmeas==1; % total ring width
   mtype='wt';
   mchar=' '; % character to go on 3rd line of header, col 62
elseif kmeas==2; % latewood width
   mtype ='wl';
   mchar='L';
elseif kmeas==3; % earlywood width
   mtype ='we';
   mchar='E';
end;



%-- Prompt for type of index
% ktype=menu('Which type of index?',...
%    'Standard index of total ring width',...
%    'Residual index of total index');
% if ktype==1; % standard index
%    dtype='std';
%    dchar=' '; % character to go on 3rd line of header
% elseif ktype==2; % residual index
%    dtype ='r';
%    dchar='R';
% end;

%-- Load the .mat input file, with data needed and the history variable
eval(['load ' pf1 ' ZI yrZI ZE yrZE Fwhen;']);

%-- Prompt for sample size variable
kss = menu('What for sample size?',...
   'Number of trees',...
   'Number of cores');
if kss==1;
   ksamp='ntree';
else;
   ksamp='ncore';
end;

%-- Prompt for 6-character site filename 
prompt={'Enter ID:'};
def={'xxxxxx'};
dlgTitle='6 Character Site ID to be used in building filename';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
id=upper(answer{1});
if length(id)~=6;
    error('id must have length 6');
end;


%-- HEADER INFORMATION APPLICABLE TO BOTH STANDARD AND RESID CHRON AND NOT DERIVED DIRECTLY FROM THE DATA

khead=menu('HEADER OPTIONS',...
    'Just simple 1-line header',...
    'ITRDB 3-line header, specify all info',...
    'ITRDB 3-line header, modify from another .crn file');
if khead==1; % simple header
    smtxt=[id,'; '];
    if strcmp(ksamp,'ntree');
        sstxt='Number of trees is used as sample size';
    elseif strcmp(ksamp,'ncore');
        sstxt='Number of cores is used as sample size';
    end;
elseif khead==2; %ITRDB 3-line header, specify all info
    hdA1=[id ' 1 ']; hdB1=[id ' 2 ']; hdC1=[id ' 3 ']; % left-hand side of each line
    
    %--- Site name. Must fit in cols 10-60.  max of 51 chars
    hdA2=blanks(51);
    prompt={'Enter site name:'};
    def={'THIS TREE-RING SITE'};
    dlgTitle='Site Name (max of 51 chars)';
    lineNo=1;
    answer=inputdlg(prompt,dlgTitle,lineNo,def);
    hdtemp=answer{1};
    ntemp = length(hdtemp);
    if ntemp>51;
        hdA2=[hdtemp(1:51) ' '];
    else;
        hdA2=[hdtemp blanks(51-ntemp) ' '];
    end;
    
    %--- Species codes
    splist={'Western Juniper',...
            'Douglas-Fir',...
            'Ponderosa Pine',...
            'Single Needle Pinyon',...
            'Sugar Pine',...
            'Limber Pine',...
            'Bigcone Douglas-fir'};
    scodes ={'JUOC','PSME','PIPO','PIMO','PILA','PIFL','PSMA'};
    kspecies=menu('Which Species',splist);
    hdA3=[scodes{kspecies} blanks(15)];
        
    %--- State (cols 10-39, 30 char max)
    hdB2=blanks(30);
    statelist={'ARIZONA','CALIFORNIA','COLORADO','MONTANA'...
            'NEVADA','NEW MEXICO','OREGON','UTAH'...
            'WASHINGTON','WYOMING',...
            'BAJA CALIFORNIA','MEXICO','CANADA'};
    kstate=menu('Where is the site?',statelist);
    hdtemp = statelist{kstate};
    ntemp=length(hdtemp);
    nfill = 30-ntemp;
    if ntemp<30;
        hdB2=[hdtemp blanks(nfill) ' '];
    else;
        hdB2=[hdtemp(1:30) ' '];
    end;
    
    %--- Elevation (meters)
    prompt={'Enter elevation:'};
    def={'0'};
    dlgTitle='Elevation of Site (meters)';
    lineNo=1;
    answer=inputdlg(prompt,dlgTitle,lineNo,def);
    zlev = str2num(answer{1});
    if zlev==0;
        kquest1 = questdlg('Is elevation really 0 m?');
        kquest1=upper(kquest1);
        if strcmp(kquest1,'YES');
            zlev=input('Elevation (m): ');
        end;
    end;
    strlev=sprintf('%4.0f',zlev);
    hdB3=[strlev 'M '];
    
    %--- Latitude/Long
    
    prompt={'Enter degrees latitude ',...
            'Enter minutes latitude',...
            'Enter N or S',...
            'Enter degrees longitude ',...
            'Enter minutes longitude',...
            'Enter W or E'};
    def={'32','00','N','110','00','W'};
    dlgTitle='Location of Tree-Ring Site';
    lineNo=1;
    answer=inputdlg(prompt,dlgTitle,lineNo,def);
    latdeg=str2num(answer{1});
    latmin=str2num(answer{2});
    latdir=answer{3};
    if strcmp(latdir,'N');
        latsign='+';
    elseif strcmp(latdir,'S');
        latsign='-';
    else;
        error('latitude must be either N or S');
    end;
    londeg=str2num(answer{4});
    lonmin=str2num(answer{5});
    londir=answer{6};
    if strcmp(londir,'W');
        lonsign='-';
    elseif strcmp(londir,'E');
        lonsign='-';
    else;
        error('latitude must be either W or E');
    end;
    
    if latdeg<0 | latmin<0 |londeg<0 | lonmin<0;
        error('Latitude and longitude units must be positive');
    end;
    if latdeg>90 | londeg>180;
        error('Latitude and longitude must be in ranges 0-90 and 0-180 degrees');
    end;
    if latmin>60 | lonmin>60;
        error('Lat and Long minutes must be in range 0-60');
    end;
    stra=sprintf('%2.0f',latdeg);
    strb=sprintf('%2.0f',latmin);
    strc=sprintf('%3.0f',londeg);
    strd=sprintf('%2.0f',lonmin);
    hdB4=[latsign stra strb lonsign strc strd blanks(9)];
    
    %--- Collector/contributor, max of 51 char, cols 10:60
    prompt={'Enter Names of Individuals:'};
    def={'MEKO, COOK, STAHLE'};
    dlgTitle='Collector/Developer/ (max of 51 chars)';
    lineNo=1;
    answer=inputdlg(prompt,dlgTitle,lineNo,def);
    hdC2=answer{1};
    ntemp=length(hdC2);
    nfill=51-ntemp;
    if ntemp>51;
        hdC2=[hdC2(1:51) ' '];
    else;
        hdC2=[hdC2 blanks(nfill) ' '];
    end;
    
elseif khead==3; % modify header off an existing .crn file
    
    [file2,path2]=uigetfile('*.crn','File to copy header lines from');
    pf2=[path2 file2];
    F = textread(pf2,'%s','delimiter','\n','whitespace','');
    F=char(F);
    F=F(1:3,:);  % pull the header lines
    
    sitename = F(1,10:60);
    speccode = F(1,62:65);
    F(3,62)=mchar;  % measurement type has been prompted for
        
    % Modify site name
    prompt={'Enter name of site:'};
    def={sitename};
    dlgTitle='Name of Site';
    lineNo=1;
    answer=inputdlg(prompt,dlgTitle,lineNo,def);
    sname = answer{1};
    ntemp=length(sname);
    nfill=51-ntemp;
    if ntemp>51;
        sname=sname(1:51);
    else;
        sname=[sname blanks(nfill)];
    end;
    F(1,10:60)=sname;
    
    % Modify species
    prompt={'Enter species code:'};
    def={speccode};
    dlgTitle='Species code (4 capital letters)';
    lineNo=1;
    answer=inputdlg(prompt,dlgTitle,lineNo,def);
    scode = upper(answer{1});
    ntemp=length(scode);
    nfill=4-ntemp;
    if ntemp>4;
        scode=scode(1:4);
    elseif ntemp==4;;
        % no action needed
    else; 
        error('Species code must be 4 letters');
    end;
    F(1,62:65)=scode;
    
end;  % if khead ==1
  
    
for nrun=1:2; % Two passes, first for standard chron, then for residual chron
    
    if nrun==1; % standard index
        dtype='std';
        dchar=' '; % character to go on 3rd line of header
        Z=ZI;
        yrZ=yrZI;
        fln = [id '.CRN']; % output filename
    elseif nrun==2; % residual index
        dtype ='r';
        dchar='R';
        Z=ZE;
        yrZ=yrZE;
        fln = [id 'r.CRN']; % output filename
    end;
     
    
    %-- Pull the tree-ring index
    X = Z(:,1);
    yrs = yrZ; % year vector for X
    
    
    % Store sample size
    if kss==1; % if using number of trees
        n = Z(:,6);
    else; % if using number of cores
        ksamp='ncore';
        n = Z(:,7);
    end;
          
    
    % Open the output file for writing
    fid=fopen(fln,'w');
    if fid==-1;
        error(['Could not open ' fln ' for writing']);
    end;
    
    % WRITE HEADERS 
    
    if khead==1; % simple header
        % Write the string data at the top of the file
        fprintf(fid,'%s   \n',[smtxt sstxt]);
    elseif khead==2; %ITRDB 3-line header, specify all info
        %--- First and last years of data
        Lgood=~isnan(X);
        yrfirst= min(yrs(Lgood));
        yrlast=max(yrs(Lgood));
        stryr1 = sprintf('%5.0f',yrfirst);
        stryr2 = sprintf('%5.0f',yrlast);
        hdB5=[stryr1 stryr2 blanks(4)];
        
        %--- Measurement and chron-type code
        hdC3=[mchar dchar blanks(17)];
        
        %-- Combine sub-headers to build lines A, B , C
        hdB=[hdB1 hdB2 hdB3 hdB4 hdB5];
        hdA=[hdA1 hdA2 hdA3];
        hdC=[hdC1 hdC2 hdC3];
        
        %--- Print the headers
        fprintf(fid,'%s\n',hdA);
        fprintf(fid,'%s\n',hdB);
        fprintf(fid,'%s\n',hdC);
        
    elseif khead==3; % modify header off an existing .crn file
        
        %--- First and last years of data -- these must be derived from the data
        Lgood=~isnan(X);
        yrfirst= min(yrs(Lgood));
        yrlast=max(yrs(Lgood));
        stryr1 = sprintf('%5.0f',yrfirst);
        stryr2 = sprintf('%5.0f',yrlast);
        hdB5=[stryr1 stryr2 blanks(4)];
        
        datacode=dchar;
        F(3,63)=dchar;
        
        %--- Change header info on first and last year to match the data in this file
        F(2,67:80)=hdB5;
        
        %--- Print the headers
        fprintf(fid,'%s\n',F(1,:));
        fprintf(fid,'%s\n',F(2,:));
        fprintf(fid,'%s\n',F(3,:));
       
    end;
    
    
    % Convert all the NaN's in X into 9990's and the sample size for those years to 0
    lmx=isnan(X);
    X(lmx)=9.990;
    n(lmx)=0;
    
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
    
    fapnd=[];
    % Append 9990's in the beginning of the data if necessary
    if rem(yrs(1),10)~=0,
        fapnd=ones(rem(yrs(1),10),1)*9990;
    end
    
    bfcount=bcount+length(fapnd);
    
    % Append 9990's at the end of the data if necessary
    if rem(yrs(length(yrs)),10)~=0,
        rapnd=ones(9-rem(yrs(length(yrs)),10),1)*9990;
    end
    rfcount=rcount+length(rapnd);
    
    % Append 9990's in the beginning of the data if necessary
    brem=rem(bfcount,10);
    if brem~=0,
        bwhln=fix(bfcount/10);
        if bwhln==0,
            fprintf(fid,'%s',id);
            fprintf(fid,'%4d',yrs(1));
            for i=1:length(fapnd),
                fprintf(fid,'%4d%3d',fapnd(i),0);
            end
        else
            fprintf(fid,'%s',id);
            fprintf(fid,'%4d',yrs(bwhln*10+brem-length(fapnd)+1));
            for i=1:brem,
                fprintf(fid,'%4d%3d',9990,0);
            end
        end
    end
    
    disp(['Hang on -- writing ' fln ]);
    
    % Main loop : write the actual data X
    for i=bcount+1:length(X)-rcount,
        if rem(yrs(i),10)==0,
            fprintf(fid,'\n%s%4d',id,yrs(i));
        end
        fprintf(fid,'%4d%3d',round(X(i)*1000),n(i));
    end
    
    % Append 9990's at the end of the data if necessary
    rrem=rem(rfcount,10);
    if rrem~=0,
        for i=1:rrem,
            fprintf(fid,'%4d%3d',9990,0);
        end
    end
    
    fclose(fid);  	% Close the input file
    
end; % for nrun  

% Update the history file to record that these are the most recent .crn files created from the .mat storage file
ctime=clock;
ctime=num2str(ctime(4:5));
dtime=date;
Fwhen{7,1}='itrdbcrn';
Fwhen{7,3}=file1;
Fwhen{7,4}=[id '.crn'];
Fwhen{7,2}=[dtime ', ' ctime];
eval(['save ' pf1 ' Fwhen -append;']);

