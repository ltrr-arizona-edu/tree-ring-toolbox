% SKEMOD.M
% This is a script .m file to modify an existing .ske file to 
% correct for false or locally absent rings. Reads data from 
% the input file and write the modified data onto the output  
% file.
% 
% USER WRITTEN FUNCTIONS NEEDED
% 
% SKPRCS.M 	Divides the input data into three blocks 
% SLVMNU.M	A modified menu function
% USINP.M	A modified menu function
% JDISP.M	A non-interactive display function
% UDISP.M	An interactive display function
% SHLMNU.M	A modified menu function
%_______________________________________________________________

global fln;
global fwn;

% Change to the default directory
cd c:\mlb\data;

% prompt for data file name. Then add suffix .dat
fln=uigetfile('*.sk*','Enter SKE filename');

% If Cancel button is pressed in the input window, flag is 
% set to -1 and control is returned to calling function
if fln==0,
  fl=-1;
  return;
end 

% Check for existence of the input file name
while ~exist(fln),
  fln=uigetfile('*.sk*','Enter Correct SKE filename');
end

lfln=length(fln);
skn=fln(1:lfln-4);
fxt=str2num(fln(lfln));
if isempty(fxt),
  fxt=1;
else
  fxt=fxt+1;
end

fwn=[skn,'.SK',num2str(fxt)];

% Open the data file for reading
fid=fopen(fln,'r');

% Read the data file in %g format and store the data 
% in a vector
sv=fscanf(fid,'%g',inf);
tm=length(sv);
m=tm/2;         % Number of rows in the data matrix

% Form an (m x n) matrix including the years 
st=reshape(sv,2,m);
st=st';
mxs=max(st(:,1));

% Plot the old data
h1=figure('Position',[50 50 550 300]);
plot(st(:,1),st(:,2),'bo');

% Testing if the file contains any zero index
ltst=st(:,2)==0;
if ~isempty(ltst(ltst)),
   udisp('WARNING : File Contains Zero Ring Index !!');
end

% Open a menu window asking the choose one out of 5 options
acmv=['Add LA Ring   ';
      'Remove LA Ring';
      'False to Real ';
      'Real to False ';
      'Quit          '];
acmn=slvmnu('Choose action',acmv);

% Ask for year index to be modified and then prompt for the 
% corresponding data values
if acmn==1,
  hsp=jdisp('Please return to the Command Window for input');
  rnb=input(['Add at ring #  ( ',num2str(st(1,1)),'-',num2str(st(m,1)),' ) = ']);

  [yrs,yrf,sts,stf,chnd] = skprcs(rnb,st,m,1);  

  cskv=input('Please enter the value corresponding to t = ');
  cskv1=input('Please enter the value corresponding to t+1 = ');
  close(hsp);
  stmm=[cskv;cskv1];
  yrmm=[rnb;rnb+1];
  lst=stmm==0;
  stmm(lst)=[];
  yrmm(lst)=[];
  dsy=[sts;stmm;stf];
  yv=[yrs;yrmm;yrf];
elseif acmn==2,
  hsp=jdisp('Please return to the Command Window for input');
  rnb=input(['Delete Ring #  ( ',num2str(st(1,1)),'-',num2str(st(m,1)),' ) = ']);

  [yrs,yrf,sts,stf,chnd] = skprcs(rnb,st,m,2);

  cskv=input('Please enter the value corresponding to t = ');
  close(hsp);
  stmm=cskv;
  yrmm=rnb;
  lst=stmm==0;
  stmm(lst)=[];
  yrmm(lst)=[];
  dsy=[sts;stmm;stf];
  yv=[yrs;yrmm;yrf];
elseif acmn==3,
  hsp=jdisp('Please return to the Command Window for input');
  rnb=input(['Real to False Ring #  ( ',num2str(st(1,1)),'-',num2str(st(m,1)),' ) = ']);
  
  [yrs,yrf,sts,stf,chnd] = skprcs(rnb,st,m,3);

  cskv=input('Please enter the value corresponding to t = ');
  cskv1=input('Please enter the value corresponding to t+1 = ');
  cskv2=input('Please enter the value corresponding to t+2 = ');
  close(hsp);
  stmm=[cskv;cskv1;cskv2];
  yrmm=[rnb;rnb+1;rnb+2];
  lst=stmm==0;
  stmm(lst)=[];
  yrmm(lst)=[];
  dsy=[sts;stmm;stf];
  yv=[yrs;yrmm;yrf];
elseif acmn==4,
  hsp=jdisp('Please return to the Command Window for input');
  rnb=input(['Real to False Ring #  ( ',num2str(st(1,1)),'-',num2str(st(m,1)),' ) = ']);

  [yrs,yrf,sts,stf,chnd] = skprcs(rnb,st,m,4);

  cskv=input('Please enter the value corresponding to t = ');
  cskv1=input('Please enter the value corresponding to t+1 = ');
  close(hsp);
  stmm=[cskv;cskv1];
  yrmm=[rnb;rnb+1];
  lst=stmm==0;
  stmm(lst)=[];
  yrmm(lst)=[];
  dsy=[sts;stmm;stf];
  yv=[yrs;yrmm;yrf];
elseif acmn==5,
  close all;
  return;
end

var=[yv dsy];

% Plot the old and new SKE data 
figure(h1);
hold on;
plot(var(:,1),var(:,2),'rx');
hold on;
plot(var(:,1),var(:,2),'r');
title(['O/X = OLD/NEW, ',fln,' / ',fwn,', ',acmv(acmn,:)]);
xlabel('Year');ylabel('Index');

chnd=min(chnd,m);
tx1=text('Position',[0.7,0.93],'String','First       Last          N',...
     'Units','normalized','Fontsize',8,'Color','b');
tx2=text('Position',[0.65,0.89],'String',['OLD      ',num2str(st(1,1)),'          ',...
    num2str(st(m,1)),'          ',num2str(m)],'Units','normalized','Fontsize',8);
tx3=text('Position',[0.65,0.85],'String',['NEW     ',num2str(yv(1)),'          ',...
   num2str(yv(length(yv))),'          ',num2str(length(yv))],'Units','normalized','Fontsize',8);

% Zooming capability
kzm=1;
vxs=axis;
while kzm<3,
  kzm=shlmnu('Choose one','Zoom in','Zoom out','Quit');

  if kzm==1,
    dh1=jdisp('Please click the mouse at the two desired corner points');
    figure(h1);
    [nxs,nys]=ginput(2);
    axis([min(nxs) max(nxs) min(nys) max(nys)]);
    grid on;
    close(dh1);
  elseif kzm==2,
    figure(h1);
    axis(vxs);
  end
end

% Open a file to put the new SKE data in
[fiw,mes]=fopen(fwn,'w');

[vr,vc]=size(var);

for i=1:vr,
  for j=1:vc,
    fprintf(fiw,'%5g',var(i,j));
  end
  fprintf(fiw,'\n');
end

fclose('all'); % Close all the open files

% Prompt for closing all graph windows
knm=usinp('Close the graph ?');
if knm,
  close all;
end 

% End of file

