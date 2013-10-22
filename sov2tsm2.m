function sov2tsm2(flnlist,whatvars,Jsers,yrswant,xyfile)
% sov2tsm2: build single tsm from selected series in sov's of multiple files 
% CALL: sov2tsm2(flnlist,whatvars,Jsers,yrswant,xyfile);
%
% Used to build matrix of tree indices as predictors in San Pedro River
% Reconstruction.  Have 30+ tree sites, each with its own .mat storage
% file.  In each file is a sov of tree indices, and identifying info.
% Want to cull selected series from each site.
%
% Meko 3-25-97
%
%******************** IN **********************************
%
% flnlist: ascii file of path\filenames  of tree sites
% whatvars: cell variable with three elements:
%   1 - name of strung-out vector (e.g.,  'ET')
%	 2 - name of matrix of first year, last year, and starting index
%		in ET of each series  (e.g., 'ETyrs')
%	 3 - name of string matrix with tree names (e.g., 'Tnms')
% Jsers (mJsers x nJsers); pointer to rows of item 2 above indicating
%		which series to cull for each site (e.g., [1 2 NaN NaN; 1 3 NaN NaN],
%		says series 1 and 2 for site 1; series 1 and 3 for site 2)
%		Right filled with NaN
%		Number of rows of Jsers should equal number of files in flnlist. 
%			This will be checked within sov2tsm2.
%      *** see jsersbld.m for quick way to initialize Jsers ******
% yrswant (1 x 2)i  first, last year of desired output tsm
% xyfile: ascii file of x,y mapping coordinates of siteslisted in flnlist
%		2 cols
%
%
%*********************** OUT **************************************
%
% .MAT file with tsm stored as X and year vector as yr
% ascii file with summary information on each series included in X:
%
%		col, file, treename yeargo  yearsp
%
% ascii file with x,y, coordinates corresponding to each tree in col of X
%
%******************  NOTE ***********************************************
%
% Calls sov2mat.m to get 
%
%*********************  FUNCTIONS CALLED *******************************
%
% sov2mat.m  -- build single tsm for a given .mat file

anotty = NaN;
clc


% How many variables will be in the master storage
numvars=0;
[mJsers,nJsers]=size(Jsers);

NV=anotty(ones(mJsers,1),:);  % will hold number of series for each of mJsers tree sites
for n = 1:mJsers;
   NV(n)=sum(~isnan(Jsers(n,:)));
   numvars=numvars+sum(~isnan(Jsers(n,:)));
end

Cxy=anotty(ones(numvars,1),ones(2,1)); % will hold x,y coordinates of trees in
%		X
% ICD = anotty(ones(numvars,1),:); % will hold site-sequence vector 
%	cross-referencing each tree to the 1,2,..nfiles sites. One use for
%	this vector is in telling spatrec1.m which local climate series to
%  use when the tree variables in X are used in reconstruction


% Initialize master storage array Zout and its year vector yrZout
yrZout=[yrswant(1):yrswant(2)]';
Zout=anotty(ones(length(yrZout),1),ones(numvars,1));

% Open file with list of filenames
fid1=fopen(flnlist,'r');

%********  Count the number of files
clc
disp('Counting the .mat files in the filelist...');
nfiles=0;
k1=1; % while control
while k1;
  	c = fgetl(fid1); 
	if ~feof(fid1);
		nfiles=nfiles+1;
	else
		if any(isletter(c));
			nfiles=nfiles+1;
		end
   	k1=0;
	end
end; % while k1
frewind(fid1);
if nfiles~=mJsers;
   error('number of .mat files not equal to row size of Jsers');
end
disp(['   Number of .mat files = ' int2str(nfiles)]);



%************* Size matrices to store path\filenames, and prefixes of .mat files
blnks25=blanks(25);
FILE2=blnks25(ones(nfiles,1),:);
PF2= blnks25(ones(nfiles,1),:);




%******************** Initial readthru of .mat files to make string matrices of file names
% and of combined path\filenames
disp('Making string matrix of path\filenames ...');

lenpf=0; % initialize max length of any path\filename
lenfn=0; % initialize max length of any filename, without .mat

for n=1:nfiles;
   % Get the .mat filename and its path 
   c=fgetl(fid1);
   pf2=strtok(c);
   len1=length(pf2);
   fslash=findstr(pf2,'\');
   if isempty(fslash); % no path prefix in file name ; assume .mat file in current dir
      file2=strtok(pf2,'.');
      pf2=[eval('cd') '\' pf2];
   else; % separate path and filename
      file2=pf2((max(fslash)+1):len1);
      file2=strtok(pf2,'.');
      %path2=pf2(1:max(fslash));
      pf2=pf2;
   end
   lenpf=max(lenpf,length(pf2));
   lenfn=max(lenfn,length(file2));
   PF2(n,1:length(pf2))=pf2;
   FILE2(n,1:length(file2))=file2;
   sprintf('%s',[int2str(n) pf2]);
end

% Trim unneeded rightmost cols off FILE2 and PF2
FILE2=FILE2(:,1:lenfn);
PF2=PF2(:,1:lenpf);

% Close filename file
fclose(fid1);


%***************************  PRELIM LOOP TO CHECK EXISTENCE OF DESIRED VARIABLES
%  IN THE MAT FILES
disp('Looping over .mat files to check that desired variables exist ...');

% pull the variable names out of cells
vbl1=whatvars{1}; % sov
vbl2=whatvars{2}; % names
vbl3=whatvars{3}; % years info

for n=1:nfiles;
   pf2=strtok(PF2(n,:));
   
      
   % Initially clear key variables from workspace
   eval(['clear ' vbl1 ' '  vbl2  ' '   vbl3]);
   
   % Load the .mat file
   eval(['load ' pf2]);
   
   % Check for existence of desired variables
   if exist(vbl1)~=1 | exist(vbl2)~=1 | exist(vbl3)~=1;
      disp(pf2);
      error('Required variable not in above file')
   end
   
   % Which cols of this file to use; and check whether inconsistent with avail num of series
   thesevar = Jsers(n,:);
   thesevar=thesevar(~isnan(thesevar));
   if (size(eval(vbl2),1)<max(thesevar)) | (size(eval(vbl2),1)<length(thesevar));
      disp(pf2);
      error('Jsers content inconsistent with above mat file');
   end
         
end


%Compute col-range in Zout for each .mat file's series 
% know already that numvars is the total number of variables, and that NV is a cv
% of number of series for each site
ncum =cumsum(NV)+1;
ncum(length(ncum))=[];
colgo=[1; ncum];
colsp=colgo+NV -1;



%**************** LOOP OVER FILES TO GRAB DATA AND BUILD OUTPUT MATRIX
disp('Looping over .mat files to grab data and build master tsm ...');

for n = 1:nfiles;
   jget=Jsers(n,:);
   jget=(jget(~isnan(jget)))'; % cv : source variables within the sov
   jput=[colgo(n):colsp(n)]; % rv of target cols in Zout
   pf2=strtok(PF2(n,:));
   eval(['load ' pf2]);% Load the .mat file
   
   % Get tsm of data for desired series, one file
   [Xtemp,yrtemp]=sov2tsm(eval(vbl1),eval(vbl3),yrswant,jget);
   
   % Put the data in correct slot
   Zout(:,jput)=Xtemp;
      
      
end

%**************  Build the cross-reference col vector ICD ***************

isp=cumsum(NV);
igo=isp-NV+1;

for n=1:nfiles;
	nrows=NV(n);
	ICD(igo(n):isp(n),:)=n(ones(nrows,1),:);
end



%********  Save tsm year vector, and site-sequence cv in .mat file as X,yr,ICD
X=Zout;
yr=yrZout;
[file4,path4]=uiputfile('*.mat','Store the tsm and year vector in this .mat file');
pf4=[path4 file4];
eval(['save ' pf4 ' X yr ICD']);


%******************** BUILD ASCII INFO FILE *******************
disp('Building the output ascii information file ...');

[file3,path3]=uiputfile('*.txt','Ascii output information table');
pf3=[path3 file3];
fid3=fopen(pf3,'w');

fprintf(fid3,'%s\n','Trees in Master Matrix');
fprintf(fid3,'%s\n\n',['Filenames in ' flnlist]);
fprintf(fid3,'%s\n\','   Col      File       Treeseq Treeid     Years');

% Loop over .mat files
jcount=1;  % sequential counter of trees, over all sites
for n=1:nfiles;
   pf2=strtok(PF2(n,:));
   file2=FILE2(n,:);
   eval(['load ' pf2]);% Load the .mat file
   str2=sprintf('%3.0f-%s\t',n,file2); % file # and file name
   
   
   % Loop over series to be used in the file
   jcol = Jsers(n,:);
   jcol=jcol(~isnan(jcol));
   nsers=length(jcol); % number of series from this file
   for m=1:nsers;
      str1=sprintf('%4.0f',jcount); % sequential tree counter, over all files
      str3=sprintf('(%3.0f)\t',jcol(m));   % seq number of tree within its mat file
      eval(['treename = ' vbl2 '(jcol(m),:);']);
      str4=sprintf('%s\t',treename);
      eval(['years = ' vbl3 '(jcol(m),1:2);']); % start and end year
      str5=sprintf('%4.0f  %4.0f',years);
      
      strall=[str1 str2 str3 str4 str5];
      fprintf(fid3,'%s\n',strall);
      jcount=jcount+1;
   
   end
end

disp('Successful completion');
fclose (fid3);


%******************* BUILD ASCII FILE OF X,Y COORDINATES OF INDIVIDUAL TREES
disp('Building ascii x,y coordinate file for trees represented by cols of X');
%     (all tree from same .mat file (derived from same .rwl file) assigned
%		same coordinates)

% load ascii file with site coords
fid5=fopen(xyfile,'r');

% Read the presumably 2-col matrix of x,y values in as a cv
xydata=fscanf(fid5,'%g');
fclose(fid5);


% Check that xydata is consistent with number of files
len1=length(xydata);
len2=nfiles*2;
if len1 ~= len2;
	error(['site xy file has ' int2str(len1) ' points, expected ' int2str(len2)]);
end

% Separate out x and y coordinates; make into 2-col x,y matrix
iy=(2:2:len2);  % rows with y coordinates
ix = iy-1; % rows with x coords
Sxy = [xydata(ix) xydata(iy)]; % site coordinates


% Build index of start row, end row  in Cxy for coordinates for each site
% Recall that NV(nfiles x 1) has the number of trees in X at each site.
isp=cumsum(NV);
igo=isp-NV+1;

for n=1:nfiles;
	nrows=NV(n);
	Cxy(igo(n):isp(n),:)=repmat(Sxy(n,:),nrows,1);
end



% Write x,y coords as ascii file
[file5,path5]=uiputfile('*.dat','Ascii x,y coords of tree indices in X');
pf5=[path5 file5];
fid5=fopen(pf5,'w');
fprintf(fid5,'%7.3f %7.3f\n',Cxy');
fclose(fid5);








   
   
