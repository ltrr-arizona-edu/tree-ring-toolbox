function recinp
%
% USAGE : recinp
%   Reads a data file in '.rec' format. The file must have the following 
%   format : 3 rows above each block data matrix contain strings. Each 
%   block data matrix is of size (m x n). Saves the full data matrix, 
%   Data ID matrix, and a text matrix in a user specified .mat file.
%
% NO INPUTS
%
% NO OUTPUTS
%
% USER WRITTEN FUNCTIONS NEEDED
%------------------------------
% UDISP.M : A user interactive graphic display function
% JDISP.M : A non-interactive graphic display function
%
% FILES PRODUCED
% --------------
% ?.MAT		A .mat file with user given name. Contains 
%		a data matrix, a data ID matrix, and a text
%		matrix explaining different columns of the 
%		data ID matrix.
%________________________________________________________
    
global fln; % Must be able to find the file in global space

% Change to the default directory
cd c:\mlb\rec;

% Prompt for data file name
fln=uigetfile('*.rec','INPUT filename');

% If Cancel button is pressed in the input window, flag is 
% set to -1 and control is returned to calling function
if fln==0,
  fl=-1;
  return;
end 

% Check for existence of the input file name
while ~exist(fln),
  fln=uigetfile('*.dat','Enter Correct filename');
end
jh=jdisp(['Reading the file: ',fln,'. Please wait ... ']); 
pause(1);

% Standard data block size of .rec files
m=28;	% Row size of each block of data
n=11;	% Column size of each block of data

% Open the data file for reading
fid=fopen(fln,'r');

% Read the first 3 lines before the data; these lines contain 
% text strings

ln1=fgetl(fid);
lnt1=length(ln1);

ln2=fgetl(fid);
lnt2=length(ln2);

ln3=fgetl(fid);
lnt3=length(ln3);

% Initialize different data and string matrices

R=[];			% Initialize data matrix

lmtx1=[];
lmtx1=[lmtx1;ln1];

lmtx2=[];
lmtx2=[lmtx2;ln2];

lmtx3=[];
lmtx3=[lmtx3;ln3];

while 1, 

  % Read rest of the data file in %g format block by block 
  % and store the data block in R matrix

  zv=fscanf(fid,'%g',[n,m]);

  if isempty(zv),	% Stop reading if end of file is reached
    break;
  end

  yrv=zv(1,:);
  zv(1,:)=[];
  zv=zv(:);
  R=[R zv];

  % Read the 3 lines before the data; Put these lines in 
  % appropriate string matrices 

  ln1=fgetl(fid);
  if length(ln1)<lnt1,
    dum=[];
    for j=1:lnt1-length(ln1),
      dum=[dum,' '];
    end
    ln1=[dum,ln1];
  end
  if length(ln1)==lnt1,
    lmtx1=[lmtx1;ln1];
  end
  lnt1=length(ln1);

  ln2=fgetl(fid);
  if length(ln2)==lnt2,
    lmtx2=[lmtx2;ln2];
  end
  lnt2=length(ln2);

  ln3=fgetl(fid);
  if length(ln3)==lnt3,
    lmtx3=[lmtx3;ln3];
  end
  lnt3=length(ln3);

end

% Take out all the alphabetic characters form the data ID 
% matrices lmtx1, lmtx2, and lmtx3

[lrx1,ltx1]=size(lmtx1);
lmtx1(lrx1,:)=[];
lmtx1(:,1:18)=[];
lmtx1(:,6:12)=[];
[lrx1,ltx1]=size(lmtx1);
lmtx1(:,ltx1-15:ltx1-8)=[];
[lrx1,ltx1]=size(lmtx1);
lmtx1(:,ltx1-1:ltx1)=[];
lm1=str2num(lmtx1);

lmtx2(:,1:14)=[];
lmtx2(:,11:26)=[];
lmtx2(:,22:39)=[];
lm2=str2num(lmtx2);

lmtx3(:,1:6)=[];
lmtx3(:,7:12)=[];
lmtx3(:,15:19)=[];
lmtx3(:,24:27)=[];
lmtx3(:,31:34)=[];
lmtx3(:,37:41)=[];
lmtx3(:,45:48)=[];
lm3=str2num(lmtx3);

% Construct the data ID matrix (Numeric)
LM=[lm1,lm2,lm3];

% A text matrix containing the column info about the data 
% ID matrix LM is made as follows :

txtm=['LM - Numeric matrix of Information';...
      '     on .dat File                 ';...
      'COL 1  : GRID POINT #             ';...
      'COL 2  : LATITUDE                 ';...
      'COL 3  : LONGITUDE                ';...
      'COL 4  : CALIB BEGIN              ';...
      'COL 5  : CALIB END                ';...
      'COL 6  : VERIF BEGIN              ';...
      'COL 7  : VERIF END                ';...
      'COL 8  : RECON BEGIN              ';...
      'COL 9  : RECON END                ';...
      'COL 10 : RSQ                      ';...
      'COL 11 : ARSQ                     ';...
      'COL 12 : ERE                      ';...
      'COL 13 : PR                       ';...
      'COL 14 : SR                       ';...
      'COL 15 : RE                       ';...
      'COL 16 : CE                       ';...
      '                                  ';...
      'R - TSM, Year in Column 1, Grid   ';...
      ' point values in remaining columns'];
 
% Make the continuous year vector and make the final data matrix
yrv=(yrv(1):yrv(m)+9)';
R=[yrv R];

% Prompt for the output file name and save data matrix R, data 
% ID matrix LM and text matrix txtm in the specified .mat file
filnam=uiputfile('*.mat','Save in MAT file ?');
eval(['save ' filnam ' R ' ' LM ' ' txtm ']);

close(jh);
jhdl=jdisp(['Successful Reading of the data file ' fln]);
pause(2);
close(jhdl);	% close the jdisp window

fclose(fid);  	% close the input file


% End of file
