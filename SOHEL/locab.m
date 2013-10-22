function [A,B]=locab(Xfn)
%
% USAGE : function [A,B]=locab(Xfn)
% Computes the number of locally absent (LA) rings in a ringwidth
% file. 
% 
% INPUTS
% ------
% Xfn		File name of the ringwidth data, .rwl assumed.
%
% OUTPUTS
% -------
% A (? x 3)	Locally absent indicator.
%   Col 1	Year
%   Col 2	Sample size (# of cores with nonmissing data)
%   Col 3	Number of cores with ring locally absent
% B (? x 3)	Subset of rows of A including only the rows 
%		with at least one locally absent value
%
% USER WRITTEN FUNCTION NEEDED 
%-----------------------------
% INVTXT.M
% JDISP.M
% PLTEXT1.M
% SHLMNU.M
% USINP.M
%______________________________________________________________

if nargin~=0,
  fln=[Xfn,'.dat'];
else
  cd c:\mlb\data;
  fln=uigetfile('*.dat','File name ?');
end

% If Cancel button is pressed in the input window, flag is 
% set to -1 and control is returned to calling function
if fln==0,
  fl=-1;
  return;
end 

% open the data file for reading
fid=fopen(fln,'r');
 
% Read the total # of columns in the file from first row
nstr=fgetl(fid);
n=str2num(nstr);
if length(nstr)>2 | isstr(n),
  fl=-1;
  disp('Wrong file type : Need the correct data file name');
  return;
end
disp(['The total Number of columns read is : ' nstr]);

% Read the first line of the file ; this line contains the 
% name of the samples
zn=fgetl(fid);
l=length(zn);          % length of first line of characters
nc=str2num(zn(1:2));   % figure the column width of data matrix
if isempty(nc),
  udisp(['Incorrect file format. Assuming a column width of 8.',...
         ' Errors may occur in the menus']);
  nc=8;
end

lnc=round(l/nc);
if l<lnc*nc,
  dum=[];
  for i=1:lnc*nc-l,
    dum=[dum,' '];
  end
  zn=[zn,dum];
elseif l>lnc*nc,
  zn=zn(1:lnc*nc);
end

% Put the sample names in a string matrix in rows and 
% indicate their sequences
zname=sprintf('%8s',zn(nc+1:2*nc));
num=1;
for i=2:lnc-1,
  znam=sprintf('%8s',zn(nc*i+1:nc*(i+1)));
  zname=[zname;znam];
  num=[num;i];
end
zname=str2mat(zname,num2str(num));

% Read rest of the data file in %g format and store the data 
% in a vector
zv=fscanf(fid,'%g',inf);
tm=length(zv);
m=tm/n;         % Number of rows in the data matrix

% Form an m x n matrix including the years 
zt=reshape(zv,n,m);
zt=zt';

% Check if the # of column read is right
B=zt(:,1);
C=ones(length(B)-1,1);
LC=C==diff(B);
if LC==1,
  fl=1;
  disp('Right number of columns read');
else
  fl=-1;
  disp('The number of columns read is wrong : Execution terminated');
  return;
end

zyr=zt(:,1);
L1=zt~=-99;
L1(:,1)=[];
nmisval=[sum(L1')]';
nmil=nmisval==0;

L2=zt==0;
L2(:,1)=[];
laval=[sum(L2')]';

nmisval(nmil)=[];
laval(nmil)=[];
zyr(nmil)=[];
A=[zyr nmisval laval];

l3=laval~=0;
B=[zyr(l3) nmisval(l3) laval(l3)];

jhdl=jdisp(['Successful Reading of the data file ' fln]);
pause(2);
close(jhdl);	% close the jdisp window

% Plotting the data
[mb,nb]=size(B);
mxs=min(mb,3);
h1=figure('Color','w');
AP=A(:,3)./A(:,2);
[APX,INAP]=sort(AP);
lapx=length(APX);
APX(1:lapx-mxs)=[];
INAP(1:lapx-mxs)=[];
%l1=APX==0;
%APX(l1)=[];
%INAP(l1)=[];
%lapx=length(APX);

subplot(211);bar(A(:,1),AP,'b');
title(['Proportion of Rings Locally Absent. Site : ',fln(1:length(fln)-4)]);
xlabel('Year');ylabel('p');
if mxs>=1,
  hold on;plot(A(INAP(mxs)),APX(mxs),'r*');
  pltext1(0.1,0.95,['* - ',num2str(A(INAP(mxs),1)),' (',num2str(APX(mxs)),')']);
end
if mxs>=2,
  hold on;plot(A(INAP(mxs-1)),APX(mxs-1),'ro');
  pltext1(0.1,0.88,['o - ',num2str(A(INAP(mxs-1),1)),' (',num2str(APX(mxs-1)),')']);
end
if mxs>=3,
  hold on;plot(A(INAP(mxs-2)),APX(mxs-2),'r+');
  pltext1(0.1,0.81,['+ - ',num2str(A(INAP(mxs-2),1)),' (',num2str(APX(mxs-2)),')']);
end

subplot(212);stairs(A(:,1),A(:,2));
title('Sample Depth');
xlabel('Year');ylabel('n');

% Zooming capability
kzm=1;
figure(h1);
subplot(211);vxs1=axis;
subplot(212);vxs2=axis;
while kzm<3,
  kzm=shlmnu('Choose one','Zoom in','Zoom out','Quit');
  if kzm==1,
    dh1=jdisp('Please click the mouse at the two desired corner points');
    pause(1);
    figure(h1);subplot(211);
    [nxs,nys]=ginput(2);
    axis([min(nxs) max(nxs) min(nys) max(nys)]);
    grid on;
    subplot(212);
    axis([min(nxs) max(nxs) vxs2(3) vxs2(4)]);
    close(dh1);
  elseif kzm==2,
    figure(h1);subplot(211);
    axis(vxs1);
    subplot(212);
    axis(vxs2);
  end
end

% Point and Get Coordinates Capability
kkn=1;
while kkn==1,
  disp('Please Click the mouse button on graph window');
  kkn=shlmnu('Get a point ?','YES','NO/Quit');
  if kkn==1,
    subplot(211);
    [gx,gy]=ginput(1);
    inda=find(A(:,1)==round(gx));
    app=A(inda,3)/A(inda,2);
    invtxt(h1,gx,app);
    if ~isempty(inda),
      subplot(212);
      invtxt(h1,A(inda,1),A(inda,2));
    end
  end
end

% Ask to close the graph window
knb=usinp('Close the graph ?');
if knb,
  close all;
end

fclose(fid);  	% close the input file

% End of file
