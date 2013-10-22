% This is a script .m file to modify an existing .ske file to 
% correct for false or locally absent rings. 

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

% open the data file for reading
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

h1=figure('Position',[50 50 550 300]);
plot(st(:,1),st(:,2),'bo');

yvo=st(:,1);

acmn=shlmnu('Choose action','Add LA Ring','Remove LA Ring',...
     'False to Real','Real to False','Quit');

if acmn==1,
  rnb=input(['Add ring #  ( ',num2str(yvo(1)),'-',num2str(yvo(m)),' ) = ']);
  chnd=find(yvo==rnb);
  chd=chnd;
  if isempty(chnd),
    while (yvo(1)-rnb)>20 | (rnb-yvo(m))>20,
      rnb=input('Ring # out of range. Add after ring #  = ');
    end
    chn=find(yvo<=rnb);
    if isempty(chn),
      chnd=1;
    else
      lchn=length(chn);
      chnd=chn(lchn);
    end
  end

  if rnb>yvo(m),
    sts=st;
    stm=[];
    stf=[]; 
    yrs=yvo(1:m,1);
    yrm=rnb;
    yrf=rnb+1;
  elseif rnb<yvo(1),
    sts=[];
    stm=[];
    stf=st; 
    yrs=rnb;
    yrm=rnb+1;
    yrf=yvo(1:m,1)+1;
  else
    if isempty(chd),
      sts=st(1:chnd,:);
      stm=[];
      stf=st(chnd+1:m,:);
      yrs=yvo(1:chnd,1);
      yrm=[rnb;rnb+1];
      yrf=yvo(chnd+1:m,1)+1;
    else
      sts=st(1:chnd-1,:);
      stm=st(chnd,:);
      stf=st(chnd+1:m,:); 
      yrs=yvo(1:chnd,1);
      yrm=yvo(chnd,1)+1;
      yrf=yvo(chnd+1:m,1)+1;
    end
  end

  cskv=input('Please enter the value corresponding to t = ');
  cskv1=input('Please enter the value corresponding to t+1 = ');
  stmm=[cskv;cskv1];
  if isempty(sts),
    dsy=[stmm;stf(:,2)];
    yv=[yrs;yrm;yrf];
  elseif isempty(stf),
    dsy=[sts(:,2);stmm];
    yv=[yrs;yrm;yrf];
  else
    dsy=[sts(:,2);stmm;stf(:,2)];
    yv=[yrs;yrm;yrf];
  end
elseif acmn==2,
  rnb=input(['Delete Ring #  ( ',num2str(yvo(1)),'-',num2str(yvo(m)),' ) = ']);
  chnd=find(yvo==rnb);
  while isempty(chnd),
    rnb=input('Incorrect #, delete ring #  = ');
    chnd=find(yvo==rnb);
  end   
  sts=st(1:chnd-1,:);
  stm=st(chnd:chnd+1,:);
  stf=st(chnd+2:m,:); 
  cskv=input('Please enter the value corresponding to t = ');
  cskv1=input('Please enter the value corresponding to t+1 = ');
  stmm=[cskv;cskv1];
  lst=stmm==zeros(2,1);
  stmm(lst)=[];
  dsy=[sts(:,2);stmm;stf(:,2)];
  ikv=0;
  if cskv==0,
    st(chnd,:)=[];
    ikv=ikv+1;
  end
  if cskv1==0,
    st(chnd+1,:)=[];
    ikv=ikv+1;
  end
  st(chnd+ikv:m-ikv,1)=st(chnd+ikv:m-ikv,1)-1;
  yv=st(:,1);
elseif acmn==3,
  rnb=input(['Real to False Ring #  ( ',num2str(yvo(1)),'-',num2str(yvo(m)),' ) = ']);
  chnd=find(yvo==rnb);
  
  close all;
  return;
elseif acmn==4,
  close all;
  return;
elseif acmn==5,
  close all;
  return;
end

var=[yv dsy];

figure(h1);
hold on;
plot(var(:,1),var(:,2),'rx');
hold on;
plot(var(:,1),var(:,2),'r');
title(['O/X = OLD/NEW, ',fln,' / ',fwn]);
xlabel('Year');ylabel('Index');

chnd=min(chnd,m);
tx1=text('Position',[0.7,0.9],'String','First       Last          N',...
     'Units','normalized','Fontsize',8,'Color','b');
tx2=text('Position',[0.65,0.86],'String',['OLD      ',num2str(yvo(1)),'          ',...
    num2str(yvo(m)),'          ',num2str(yvo(chnd))],'Units','normalized','Fontsize',8);
tx3=text('Position',[0.65,0.82],'String',['NEW     ',num2str(yv(1)),'          ',...
   num2str(yv(length(yv))),'          ',num2str(yv(chnd+1))],'Units','normalized','Fontsize',8);

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

[fiw,mes]=fopen(fwn,'w');

[vr,vc]=size(var);

for i=1:vr,
  for j=1:vc,
    fprintf(fiw,'%5g',var(i,j));
  end
  fprintf(fiw,'\n');
end

fclose('all');

knm=usinp('Close the graph ?');
if knm,
  close all;
end 

% End of file

