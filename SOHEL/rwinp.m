function [zname,zyrs,z,zind,fl]=rwinp
%
% USAGE : [zname,zy,z,zind,fl]=rwinp
%  Reads a data file in '.dat' format:
%
%	Row 1 contains series names (max of 8 characters, starting
%	with a letter, separaterd by blanks
%
%	Remaining rows have the data -- a year col and one col per 
%	series
%
%	For compatibility with PROGRAM LIBRARY "YUX" output, row 1
%	might have an initial string beginning with a number, like
%	"3_cols".  This string is skipped by rwinp.m
%
% NO INPUTS
%
% OUTPUTS
%--------
% zname		An string matrix of names (8cols) and sequence (3 cols) 
% zyrs (? x 3)	A 3-column matrix of first and last year for each series. 
% z (? x 1)	The data array in a strung-out column vector.
% zind (? x 1)	The starting indices of each series in z in a column vector.
% fl (1 x 1)	fl = 1 if file is successfully read.
%		fl = -1 if file could not correctly read.
%
% USER WRITTEN FUNCTIONS NEEDED
%------------------------------
% UDISP.M : A user interactive graphic display function
% JDISP.M : A non-interactive graphic display function
%________________________________________________________
    
global fln; % Must be able to find the file in global space

% Change to the default directory
cd c:\mlb\data;

% prompt for data file name
fln=uigetfile('*.dat;*.sat','INPUT filename');

% If Cancel button is pressed in the input window, flag is 
% set to -1 and control is returned to calling function
if fln==0,
  fl=-1;
  return;
end 

% Check for existence of the input file name
while ~exist(fln),
  fln=uigetfile('*.dat','Enter Correct RW filename');
end

% open the data file for reading
fid=fopen(fln,'r');
 



% Read the first line of the file ; this line contains the 
% names of the series
% zname will have series names in first 8 cols, seq nos in
%	cols 9-11
zn=fgetl(fid);
[zname,nsers]=namenum(zn);
disp(['The number of series is : ',int2str(nsers)]);
n=nsers+1;  % a years column added

% Read rest of the data file in %g format and store the data 
% in a vector
zv=fscanf(fid,'%g',inf);
tm=length(zv);
m=tm/(n);         % Number of rows in the data matrix

% Form an m x n matrix including the years 
zt=reshape(zv,n,m);
zt=zt';

% Check if the # of column read is right
B=zt(:,1);
C=ones(length(B)-1,1);
LC=C==diff(B);
if all(LC),
  fl=1;
  disp('Year col increments by 1, as it should');
else
  fl=-1;
  disp('Error: year col does not increment by 1');
  return;
end

% Initialize the return variables
zyr=zeros(m,1);
zyrs=zeros(n-1,3);
z=zeros(m*(n-1),1);
L=z;
zind=zeros(n,1);
nzisum=1;

zyr=zt(:,1);
L1= (zt~=-99) & (~isnan(zt)) ;
L1(:,1)=[];

zind(1)=nzisum;
for i=1:n-1,
  % Put the beginning and ending years in a matrix
  minyr = min(zyr(L1(:,i)));
  maxyr = max(zyr(L1(:,i)));
  zyrs(i,:)=[i minyr  maxyr];
  % Put the starting indices of each series in a column vector
  nzisum = nzisum+ (maxyr-minyr+1);
  zind(i+1)=nzisum;
end

z(:)=zt(:,2:n);  % Put the modified columns in a column vector z
L(:)=L1;
z(~L)=[];         % Delete all -99 values from each column

jhdl=jdisp(['Successful Reading of the data file ' fln]);
pause(2);
close(jhdl);	% close the jdisp window

kk=[];
kk=usinp('View the series info ?');
if isempty(kk),
  kk=0;
end
if kk~=0,
  [ml,nl]=size(zname);
  [m,n]=size(zyrs);
  zyr=[];
  znamm=zname; % Do not change zname matrix
  znamm(:,nl)=[];

  for i=1:m,
    zyrs1=num2str(zyrs(i,1));
    zyrs2=num2str(zyrs(i,2));
    zyrsn=num2str(zyrs(i,n));
    lzy=[length(zyrs1),length(zyrs2),length(zyrsn)];
    if (lzy(1)|lzy(2)|lzy(3))<4 & max(lzy)<5,
      for j=1:(4-lzy(1)),
        zyrs1=[zyrs1,' '];
      end
      for j=1:(4-lzy(2)),
        zyrs2=[zyrs2,' '];
      end
      for j=1:(4-lzy(3)),
        zyrsn=[zyrsn,' '];
      end
    end
    zyr=[zyr;zyrs1,'   ',znamm(i,:),'   ',...
         zyrs2,'   ',zyrsn];
  end

  disp(zyr)    % View the series information.
  %disp('Press any key to continue');
  pause(2);
end	% end of if loop

fclose(fid);  	% close the input file

% End of file
