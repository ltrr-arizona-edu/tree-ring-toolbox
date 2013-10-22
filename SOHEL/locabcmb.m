%function [A,B]=locabcml(filnam)
%
% USAGE : Returns A and B matrices containing cumulative 
%         nonmissing and locally absent tree ring data from 
%	  different sites
%
% INPUTS 
%-------
% filnam	A text string containing the ascii file name
%		without the extension .cat.
% OUTPUTS
% -------
% A (? x 3)	Cumulative Locally absent indicator.
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

%if nargin~=0,
%  fln=[filnam,'.cat'];
%else
  cd c:\mlb\data;
  fln=uigetfile('*.cat','File name ?');
  % If Cancel button is pressed in the input window, flag is 
  % set to -1 and control is returned to calling function
  if fln==0,
    fl=-1;
    return;
  end 
%end

% open the data file for reading
fid=fopen(fln,'r');

namat=[];
while 1
  line = fgets(fid);
  if ~isstr(line), break, end
  namat=[namat;line];
end
disp(namat)
fclose(fid);

[na,nb]=size(namat);

rnl=zeros(na,1);
yrbeg=zeros(na,1);
yrend=zeros(na,1);
for i=1:na,
  rnl(i)=length(namat(i,:));
  yrbeg(i)=str2num(namat(i,rnl(i)-13:rnl(i)-9));
  yrend(i)=str2num(namat(i,rnl(i)-4:rnl(i)));
end

ybg=min(yrbeg);
ynd=max(yrend);
A=zeros(ynd-ybg+1,3);
B=zeros(ynd-ybg+1,3);
A(:,1)=(ybg:ynd)';
B(:,1)=(ybg:ynd)';

for i=1:na,
  ofnam=namat(i,1:rnl(i)-15);
  L1=ofnam'==' ';
  ofnam(L1)=[];
  [AT,BT]=locabn(ofnam);
  lnat=length(AT(:,1));

  yr1=AT(1,1);
  yr2=AT(lnat,1);
  lga=A(:,1)>=yr1 & A(:,1)<=yr2;
  A(lga,2:3)=A(lga,2:3)+AT(:,2:3);

end

la=find(A(:,2)==0);
A(la,:)=[];
B=A;
lb=find(B(:,3)==0);
B(lb,:)=[];

% Plotting the data
[mb,nb]=size(B);
mxs=min(mb,3);
h1=figure('Color','w');
AP=A(:,3)./A(:,2);
[APX,INAP]=sort(AP);
lapx=length(APX);
APX(1:lapx-mxs)=[];
INAP(1:lapx-mxs)=[];

subplot(211);bar(A(:,1),AP,'b');
title(['Proportion of Rings Locally Absent. Group File : ',fln(1:length(fln)-4)]);
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
    invtxt(h1,gx,gy);
    inda=find(A(:,1)==round(gx));
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


% End of file



