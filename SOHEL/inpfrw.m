%
% USAGE : SKEINP
%   Script file skeinp.m reads data from .SKE files and 
%   writes the data in .DAT format to an output file.
%
%
% USER WRITTEN FUNCTIONS NEEDED
%------------------------------
% JDISP.M	A non-interactive window display function
%________________________________________________________

global fln;
global fwn; 

% Prompt for the Number of files to be read
hk=jdisp(['Please Return to Command Window to enter the ',...
          '# of files to be read.']);
nf=input('Number of files to be read =  ');
close(hk);

% Change to the default directory
cd c:\mlb\dave;

mxs=zeros(nf,1);

for i=1:nf,
  % prompt for data file name. Then add suffix .dat
  fln=uigetfile('*.ske','Enter SKE filename');

  % If Cancel button is pressed in the input window, flag is 
  % set to -1 and control is returned to calling function
  if fln==0,
    fl=-1;
    return;
  end 

  % Check for existence of the input file name
  while ~exist(fln),
    fln=uigetfile('*.ske','Enter Correct SKE filename');
  end
  skn=fln(1:length(fln)-4);
  eval(['skn',num2str(i),'=skn;']);

  % open the data file for reading
  fid=fopen(fln,'r');
  sknm1=fgetl(fid);
  sknm2=fgets(fid);
  iskm=max(length(sknm1),length(sknm2));
  iskn=min(length(sknm1),length(sknm2));
  skpc=' ';
  for k=1:iskm-iskn,
    skpc=[skpc,' '];
  end  

  if length(sknm1)<length(sknm2),
    sknm=[sknm1,skpc;sknm2,' '];
  elseif length(sknm1)>length(sknm2),
    sknm=[sknm1,' ';sknm2,skpc'];
  else
    sknm=[sknm1,' ';sknm2,' '];
  end

  % Read the data file in %g format and store the data 
  % in a vector
  sv=fscanf(fid,'%g',inf);
  tm=length(sv);
  m=tm;         % Number of rows in the data matrix

  % Form an m x n matrix including the years 
  st=reshape(sv,1,m);
  st=st';
  mxs(i)=max(st);
  eval(['st',num2str(i),'=st;']);
end

ms=max(mxs);

dsy=[];

for i=1:nf,
  sy=eval(['st',num2str(i)]);
  syr=length(sy);
  if syr<ms,
    sy=[sy;-99*ones(ms-syr,1)];
  end
  dsy=[dsy sy];
end
%   eval(['dsy',num2str(i),'=dsy;']);

eval([skn(1:3),'=dsy;']);

fwn=uiputfile('*.rwn','Output File Name');

fiw=fopen(fwn,'w');

[vr,vc]=size(dsy);
sknm=[sknm,'  '];

fprintf(fiw,'%d\n',vc);
fprintf(fiw,'%s\n',sknm);
for i=1:vr,
  for j=1:vc,
    fprintf(fiw,'%5g',dsy(i,j));
  end
  fprintf(fiw,'\n');
end

% End of file