function Z=canatmp1(pf1)
% canapcp1:  monthly canadian temperature series to .mat file
% Z=canatmp1(pf1);
% Last revised 4/14/00
%
% Converts a Canadian monthly temperature series from the raw format as accessed over the
% web to a .mat file
%
%*** INPUT 
%
% pf1 (1 x ?)s  path/filename to the input station monthly T series
%
%
%*** OUTPUT
%
% Z (? x 13)r  monthly temperature, degrees C
%
%
%*** REFERENCES -- none
%*** UW FUNCTIONS CALLED -- none
%*** TOOLBOXES NEEDED -- none
%
%*** NOTES
%
% Refer to the README.TXT file from the web for details on properties of the data


% Build index to numeric column data
i=[1:4 6:12 16:22 26:32 36:42 46:52 56:62 66:72 76:82  86:92 96:102 106:112 116:122];

% Read the data file
file = textread(pf1,'%s','delimiter','\n','whitespace','');


% Convert from cell to string matrix
file=char(file);

% Strip header lines
file(1:2,:)=[];


% Build a blank column
nsize=size(file,1);
bc = repmat(blanks(1),nsize,1);

% Pull numeric data, and subdivide with spaces
X=file(:,i);

% Add spaces
Y1=[X(:,1:4) bc X(:,5:11) bc X(:,12:18) bc X(:,19:25) bc X(:,26:32) bc X(:,33:39) bc];
Y2=[X(:,40:46) bc X(:,47:53) bc X(:,54:60) bc X(:,61:67) bc];
Y3=[X(:,68:74) bc X(:,75:81) bc X(:,82:88)];
Y=[Y1 Y2 Y3];

% Character to numeric
Z=str2num(Y);

% Missing values to NaN
W=Z(:,2:13);
L1=W==-9999.9;
if any(any(L1));
   W(L1)=NaN;
end;
Z=[Z(:,1) W];
yr = Z(:,1);
Zmeat=Z(:,2:13);

% If any whole "missing years", insert NaN rows
yrgo = min(Z(:,1));
yrsp = max(Z(:,1));
nfull = yrsp-yrgo+1;
Z1=repmat(NaN,nfull,13);
yr1 =(yrgo:yrsp)';
i1 = yr-yrgo+1;
Z1(i1,2:13) = Zmeat;
Z1(:,1)=yr1;
Z=Z1;