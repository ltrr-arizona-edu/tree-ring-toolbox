function xytree1
% xytree1.m  -- using .dat files :
%	(1) xy coordinates for "all" (say, 668, for set A) tree ring sites
%  (2) starting and ending years for all sites
%
%  and user-prompted info on desired minimum time coverage,
%
% Make an x,y plotting file to be a post file for mapping site 
% locations.
%

% Help in file selection by specifying letter code for tree-ring set
char=input('Letter code for tree-ring matrix: ','s');

%************ GET THE XY MAPPING COORDINATES
fn5=['crnxy' char '.dat'];
txt5='crnxy?.dat? (in D:\wrk6)-- tree-ring x,y coords';
[file5,path5]=uigetfile(fn5,txt5);
eval(['load ',path5,file5]);  
eval(['fxy = '  'CRNXY' upper(char) ';']);


%**********   GET THE START AND END YEARS
fn5=['gosp' char '.dat'];
txt5='gosp?.dat? (in D:\wrk6)-- tree-ring start,end years';
[file5,path5]=uigetfile(fn5,txt5);
eval(['load ',path5,file5]);  
eval(['IYRS = ' 'GOSP' upper(char) ';']);
% Subract 8000 from values in any row where ending year >3000
L7=IYRS(:,2)>3000;
if ~isempty(L7)
	temp=IYRS(L7,:);
	temp=temp-8000;
	IYRS(L7,:)= temp;
end


%************* CHECK THAT fxy and IYRS same row size
disp(['input xy-file has ' int2str(size(fxy,1)) ' points']);
pause (1)
if size(fxy,1) ~= size(IYRS,1),
	error('xy-file and years file not same row size')
end


%************** MAKE LOGICAL VECTOR FOR SELECTED ROWS
yr1=input('Series must start this year or earlier: ');
yr2=input('Series must have good data thru this year or later');

L=IYRS(:,1)<=yr1 & IYRS(:,2)>=yr2;

% Pull the selected coordinates
Z = fxy(L,:);  % coords of selected series
disp([int2str(size(Z,1)) ' points selected']);
pause (1)

%******** Build x,y (long-lat) plot file, with EV  for use in surfer
fn5=['xy' char '*.dat'];
txt5='post file; ex: xya1.dat,xya2.dat';
[file5,path5]=uiputfile(fn5,txt5);
fid1=fopen([path5 file5],'w')
fprintf(fid1,'%7.2f%6.2f\n',Z');
fclose (fid1)  






