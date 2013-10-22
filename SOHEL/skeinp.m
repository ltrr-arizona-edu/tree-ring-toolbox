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
cd c:\mlb\data;

mxs=zeros(nf,1);
sknm=['8 -col'];

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
  sknm=[sknm,'  ',skn];
  eval(['skn',num2str(i),'=skn;']);

  % open the data file for reading
  fid=fopen(fln,'r');

  % Read the data file in %g format and store the data 
  % in a vector
  sv=fscanf(fid,'%g',inf);
  tm=length(sv);
  m=tm/2;         % Number of rows in the data matrix

  % Form an m x n matrix including the years 
  st=reshape(sv,2,m);
  st=st';
  mxs(i)=max(st(:,1));
  eval(['st',num2str(i),'=st;']);
end

ms=max(mxs);

yv=(1:ms)';
var=yv;

for i=1:nf,
  dsy=zeros(ms,1);
  sy=eval(['st',num2str(i)]);
  syr=length(sy(:,1));
  ind=zeros(syr,1);
  if syr~=ms,
    for j=1:syr,
      ind(j)=find(yv==sy(j,1));
    end
    if ind(1)~=1,
      dsy(1:ind(1)-1)=-99*ones(ind(1)-1,1);
    end
    lind=length(ind);
    if ind(lind)~=ms,
      dsy(ind(lind)+1:ms)=-99*ones(ms-ind(lind),1);
    end
    for j=1:syr,
      dsy(ind(j))=sy(j,2);
    end
  else
    dsy=sy;
  end
  var=[var dsy];
end
%   eval(['dsy',num2str(i),'=dsy;']);

eval([skn(1:3),'=var;']);

fwn=uiputfile('*.sat','Output File Name');

fiw=fopen(fwn,'w');

[vr,vc]=size(var);
sknm=[sknm,'  '];

fprintf(fiw,'%d\n',vc);
fprintf(fiw,'%s\n',sknm);
for i=1:vr,
  for j=1:vc,
    fprintf(fiw,'%5g',var(i,j));
  end
  fprintf(fiw,'\n');
end

% End of file