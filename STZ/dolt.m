F = textread('xxxxxx.crn','%s','delimiter','\n','whitespace','');
F=char(f);
F=f(1:3,:);  % pull the header lines

sitename = f(1,10:60);
speccode = f(1,62:65);
datacode=f(3,63);
yrgo=f(2,67:71);
yrsp=f(2,72:76);

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


% Modify data code
prompt={'Enter data code:'};
def={datacode};
dlgTitle='Data code (1-letter)';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
dcode = upper(answer{1});
dcode=strjust(dcode,'left');
ntemp=length(dcode);
nfill=1-ntemp;
if ntemp>1;
   dcode=dcode(1);
elseif ntemp==1;;
   % no action needed
else; 
   error('Data code must be 1 letter');
end;
F(3,63)=dcode;
