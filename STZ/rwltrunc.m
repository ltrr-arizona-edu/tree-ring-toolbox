function rwltrunc
% rwltrunc:  truncate .rwl file to exclude all lines with data earlier than specified year
% rwltrunc;
% Last revised 1-25-01
%
% Sometimes you may want to exclude earliest ringwidth data from an .rwl.  For example, the .rwl may be
% accompanying a chronology to ITRDB, the earliest data may have been excluded from the chronology, and the 
% eearliest data would only confuse users of the chronology.  This is especially true if the earliest data
% were of questionable dating quality.  rwltrunc can be run after rwlinp to lop off any part of any ring-width 
% series that contains data before a specified year.
%
%
%*** INPUT
%
% No arguments
% User prompted for earliest acceptable start year, which must be evenly divisible by 10.
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

[file1,path1]=uigetfile('*.rwl','Infile to be truncated');
pf1=[path1 file1];

[file2,path2]=uiputfile('*.rwl','Outfile to hold truncated data');
pf2=[path2 file2];



D = textread(pf1,'%s','delimiter','\n','whitespace','');

% Start year
prompt={'Enter year'};
def={'0'};
dlgTitle='Exclude all rows of data before this year';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
yrgo = str2num(answer{1});
if mod(yrgo,10)~=0;
    error('Year must be divisible by 10');
end;

E=char(D);
cyr=E(:,8:12); % should hold year
yr = str2num(cyr);

L=yr<yrgo;
if any(L);
    fid1=fopen(pf2,'w');
    
    E(L,:)=[];
    [mE,nE]=size(E);
    for n = 1:mE;
        s=E(n,:);
        fprintf(fid1,'%s\n',s);
    end;
    fclose (fid1);
else;
    h=msgbox('No data needs to be lopped off');
end;


