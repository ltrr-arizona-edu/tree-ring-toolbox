function rwlnms
% rwlnms:  change core names in rwl file and write revised data to new rwl
% rwlnms;
% Last revised 4-12-01
%
% You have a .rwl file with irregular core names.  For example, RIN09AXL, or RIN9A. 
% You want the three-letter code, followed by tree number (with 0 as first digit if 
% two-digit), followed by core letter.   rwlnms is a quick way to get a new .rwl with revise
% names.  To use rwlnms, first use DPL\SUR on the original rwl.  Then vedit the sur output
% so that the orig names are in, say, cols 1-8, and the revised names in cols 10-17, with
% space in between.  Or you might have orig names in 1-9, new in 1-19.
%
%*** INPUT 
%
% No args
% User prompted to click on files: input .rwl, .txt file with series names, output .rwl file
%
%
%*** OUTPUT
%
% .rwl file with revised names
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES


[file1,path1]=uigetfile('*.rwl','Input .rwl file');
pf1=[path1 file1];

[file3,path3]=uigetfile('*.txt','.txt file with old and revised series names');
pf3=[path3 file3];

[file2,path2]=uiputfile('*.rwl','Output .rwl file');
pf2=[path2 file2];

prompt={'Enter start and end cols for old names:','Enter start and end cols for new names:'};
def={'[1 8]','[10 17]'};
dlgTitle=['Organization of names in '  pf3];;
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
i1 = str2num(answer{1});
i2= str2num(answer{2});

% Read input .rwl file
R1 = textread(pf1,'%s','delimiter','\n','whitespace','');
R1=char(R1);
nlines=size(R1,1);
ncols=size(R1,2);
R2=R1;

% Read names file
T = textread(pf3,'%s','delimiter','\n','whitespace','');
T1=T{1};
T=char(T);
nser=size(T,1);

    
Told = T(:,i1(1):i1(2));
Tnew = T(:,i2(1):i2(2));


for n = 1: nser;
    told=Told(n,:);
    tnew=Tnew(n,:);
    for m = 1:nlines;
        d = R1(m,:);
        if strcmp(told,d(i1(1):i1(2)));
            e = strrep(d,told,tnew);
            R2(m,:)=e;
        end;
    end;
end;

fid2=fopen(pf2,'w');
for m =1:nlines;
    s=R2(m,:);
    fprintf(fid2,'%s\n',s);
end;

fclose(fid2);



disp('here');