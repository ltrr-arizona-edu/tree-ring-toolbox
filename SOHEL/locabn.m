function [A,B]=locabn(Xfn)
%
% USAGE : function [A,B]=locabn(Xfn)
% Computes the number of locally absent (LA) rings in a ringwidth
% file. 
% 
% CALLING FUNCTION : LOCABCML.M
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
  fln=Xfn;
else
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
  jh=jdisp(['Incorrect file format. Assuming a column width of 8.',...
         ' Errors may occur in the menus']);
  pause(1);
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
num=num2str(1);
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
close(jh);

fclose(fid);  	% close the input file

% End of file
