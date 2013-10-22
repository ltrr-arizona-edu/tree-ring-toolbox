function ellook
%
% Graphical look at earlywood and latewood measurements
% Plots the following series
%  1. EWW and LWW 
%  2. RW  and EWW
%  3. RW  and LWW
%  4. RW  and (EWW+LWW)
%  5. EW % of RW
%
% D Meko 10-14-96
%
% No input or output args.  User is prompted for info and data files
%
%************** FILES THAT MUST BE AVAILABLE
%
% root*.txt -- info on  path to .eww, path to .rw files
%	Example:
%
%		d:\treedata\dcy\eww\
%		d:\treedata\dcy\rw\
%
% Early and late-growth files in specified path
% Total ring-width file in specified path
%
%
%************ USER WRITTEN FUNCTIONS CALLED
%
% rwread1.m

close all

% Get the paths to .eww, .lww and .rw files
[file1,path1]=uigetfile('root*.txt','file info');
pf1=[path1,file1];
fid1=fopen(pf1,'r')
% Extract the info on file names and root
epath=strtok(fgetl(fid1));
wpath=strtok(fgetl(fid1));
fclose(fid1);

% Get the root of file name for this core; get it from 
% directory with the eww files
[file2,path2]=uigetfile([epath '*.eww'],'Core early growth file');
file2=lower(file2);
nf2=length(file2);
cfile2=file2((nf2-2):nf2);
if ~all(cfile2=='eww')
	error('The file you picked is not an .eww file')
end
root=strtok(file2,'.');
% Check for an 'xe' on end of root, and remove it if it is there
nc = length(root);
if  all(root((nc-1):nc)=='xe');
	kxe=1;
	root=root(1:(nc-2));
else
	kxe=0;
end


% Get early, late and totalwidth data
E=rwread2(1,root,epath,wpath,kxe);
L=rwread2(2,root,epath,wpath,kxe);
T=rwread2(3,root,epath,wpath,kxe);

% Display first and last years of early, late and total width files
P=[E(1,1)  E(size(E,1),1) ...
	L(1,1)  L(size(L,1),1) ...
	T(1,1)  E(size(T,1),1)];
fmt1='%5.0f %5.0f\n';
S=sprintf(fmt1,P');

LL1=E(1,1)==L(1,1) & L(1,1)==T(1,1);
LL2=E(size(E,1),1)==L(size(L,1),1) & L(size(L,1),1)==T(size(T,1),1);
if ~all(LL1==LL2)
	disp(S)
	pause(4)
	error('EWW,LWW,TRW not identical periods -- see above');
end


% EWW and  LWW plot
figure(1)
plot(E(:,1),E(:,2),'y',L(:,1),L(:,2),'m');
legend('EWW','LWW')
title(['Core ID = ',root]);
ylabel('mm x 100')
xlabel('Year')

% EWW and RW
figure (2)
plot(E(:,1),E(:,2),'y',T(:,1),T(:,2),'m');
legend('EWW','RW')
title(['Core ID = ',root]);
ylabel('mm x 100')
xlabel('Year')

% LWW and RW
figure(3)
plot(L(:,1),L(:,2),'y',T(:,1),T(:,2),'m');
legend('LWW','RW')
title(['Core ID = ',root]);
ylabel('mm x 100')
xlabel('Year')

% RW and EWW+LWW
figure(4)
plot(T(:,1),T(:,2),'y',E(:,1),(E(:,2)+L(:,2)),'m');
legend('RW','EWW+LWW')
title(['Core ID = ',root]);
ylabel('mm x 100')
xlabel('Year')
zoom xon

% EWW as percentage of total RW
figure(5)
% Compute mean ratio, ignoring NaNs
R=(E(:,2) ./ T(:,2));
L1=~isnan(R);
meanrat=mean(R(L1));
plot(E(:,1),E(:,2) ./ T(:,2),...
	[T(1,1) T(size(T,1),1)],[meanrat meanrat]);
legend('EWW/RW','Mean')
title(['Core ID = ',root]);
ylabel('Ratio')
xlabel('Year')

% Test for trend in EWW/RW 
nsize=length(R(L1));
X=[ones(nsize,1) E(L1,1)];
y=R(L1);
alpha=0.01;
figure(6)
[B,BINT,H,RINT,STATS] = REGRESS(y,X,alpha);
yhat = X*B;
figure(7)
plot(E(L1,1),R(L1,1),'o',E(L1,1),yhat)
title([root ', Linear Trend in Ratio EWW/TRW'])
xlabel('Year')
ylabel('Ratio')
xa=[min(E(L1,1)) max(E(L1,1))];
ya=get(gca,'YLim');
x1=xa(1)+0.01*(xa(2)-xa(1));
y1=ya(2)-0.1*(ya(2)-ya(1));
x2=x1;
y2=ya(2)-0.15*(ya(2)-ya(1));
x3=x1;
y3=y2-0.05*(ya(2)-ya(1));
x4=x1;
y4=y3-0.05*(ya(2)-ya(1));
str1=sprintf('%8.5f',B(2)*100);
str2=sprintf('%6.4f-%6.4f',BINT(2,:)*100);
str3=sprintf('%5.3f',STATS(1));
str4=sprintf('%6.2f (%8.6f)',STATS(2:3));
str1=['Slope: ' str1 ' (per 100 years)'];
str2=['Conf. Limits (alpha=0.01):  ' str2];
str3=['R-squared: ' str3];
str4=['Overall-F (p-value): ' str4];
text(x1,y1,str1);
text(x2,y2,str2);
text(x3,y3,str3);
text(x4,y4,str4);
STATS

