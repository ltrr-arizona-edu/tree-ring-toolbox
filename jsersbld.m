function Jsers=jsersbld
% jsersbld:  build pointer matrix Jsers to be used by sov2tsm2.m
% CALL:  Jsers=jsersbld;
%
% Meko 4-16-97
%
%
%****************** IN 
%
% User points to file holding list of .mat files for tree sites.
% User responds to prompt for tree-ring variable type (e.g., ET)

clc
a=NaN;

% Get name of file holding file list
[file1,path1]=uigetfile('*.txt','Infile of list of .mat files');
pf1=[path1 file1];
fid1=fopen(pf1,'r');

%********  Count the number of files
disp(['Counting the files listed in: ' pf1]);
nfiles=0;
k1=1; % while control
while k1;
  	c = fgetl(fid1); 
   if ~feof(fid1);
      nfiles=nfiles+1;
   else
      if any(isletter(c));
         nfiles=nfiles+1;
      else
      end
      k1=0;
	end
end; % while k1
disp(['   Number of files = ' int2str(nfiles)]);
frewind(fid1);

%*************  Check that can open each file
disp('Checking that each .mat file can be opened');
for n=1:nfiles
   pf2=fgetl(fid1);
   strerr=['error([''No .mat file by the following name:'' pf2])'];
   string1=['load ' pf2];
   eval(string1,strerr);
end
disp('   No trouble opening any of the listed files');
frewind(fid1);      

%************  Find out which variable is to be worked with
prompt={'Enter the name of the tree-ring variable of interest (.e.g., ET)'};
title=' ';
lineNo=1;
def={'ET'};
answer=inputdlg(prompt,title,lineNo,def);
varnm1= char(answer);

%************ Which string matrix identifies the series of that variable
switch varnm1
case {'ET','IT'};
   id='Tnms';
   xmask='tmask';
case {'I','E'};
   id='nms';
   xmask='cmask';
otherwise
   id=[];
   xmask=[];
   error('Only acceptable responses are ET,IT, I or E');
end


%********  Build the matrix  Jsers
%
% Assume maximum average of 200 possible series per file
m1=200; % maximum possible col size of Jsers
Jsers = a(ones(nfiles,1),ones(m1,1)); % initialize Jsers as matrix of NaN

%*************  Loop over files, finding row size of id each time, and
% build matrix Jsers
disp('Filling in rows of Jsers ...');
for n=1:nfiles
   eval(['clear ' id]);
   pf2=fgetl(fid1);
   eval(['load ' pf2]);
   if exist(id)~=1 | exist(xmask)~=1;
      disp(['File ' pf2]);
      error(['No variable called ' id  ' or ' xmask 'in above .mat file']);
   else
      eval(['v= find(' xmask ');']);
      Jsers(n,1:length(v))=v';
   end
   
end


fclose(fid1);


% Strip off columns of all-NaN from Jsers
L=isnan(Jsers);
Jsers(:,all(L))=[];
 
