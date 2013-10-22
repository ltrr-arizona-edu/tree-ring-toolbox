function v2conv3
% v2conv3:  combine GHCN version 3 monthly mean max and mean min T files into monthly mean
% v2con3;
% Last revised 5-22-00
%
% Build .mat files of monthly mean temperature from GHCN .mat files of monthly mean
% minimum and maximum.  
%
%*** IN
%
% No args.
% Prompted to choose directories with the maxT files ('D'), the  
%   the minT files ('E') and the directory for the mean-mean files ('M')
%
%*** OUT
%
% No args.
% .mat files whose filenames start with 'M' are generated and stored in the 
%   specified directory
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED --NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% You have built the input maxT and minT files using v2con1.m.
% Program checks that
%  -every minT file has a matching maxT
%  -year coverage is the same
%  -missing data slots are identical
% Program then takes mean of min and max to get monthly mean
% Takes 3 min to do 1002 files

t0=clock;

prompt={'For Max T:','For Min T:'};
def={'c:\data\ghcn\tmpv2_0\matflsd\','c:\data\ghcn\tmpv2_0\matflse\'};
dlgTitle='Enter Input Directories';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
pfd=answer{1};
pfe=answer{2};


prompt={'Output directory:'};
def={'c:\data\ghcn\tmpv2_0\matflsm\'};
dlgTitle='Enter Output Directory for Computed Mean';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
pfm=answer{1};

sufx='mat';

% maxT directory and build a file list
Sd=dirfls(pfd,sufx,2); % cell with filename strings

% minT directory and build a file list
Se=dirfls(pfe,sufx,2); % cell with filename strings

% Check that same number of .mat files in the max and min dirs
[m1,n1]=size(Sd);
[m2,n2]=size(Se);
if m1~=m2 | n1~=n2;
   error([pfd ' and ' pfe ' have different numbers of .mat files']);
end;

% Check that all maxT files start with d, all min with e
if ~all(Sd(:,1)=='d');
   error('Not all maxT files start with d');
end;
if ~all(Se(:,1)=='e');
   error('Not all minT files start with e');
end;

% Check that file names of maxT and minT files are identical except for first letter
s1=Sd(:,2:n1); % lops off leading letter
s2=Se(:,2:n2);
% Convert to cell arrays of strings, then sort
s1a=sort(cellstr(s1)); s2a=sort(cellstr(s2));
% Convert sorted cells back to string matrix
s1b=char(s1a); s2b=char(s2a);
if any(any(s1b~=s2b));
   error('Filenames (except for first letter) of minT and maxT not identical');
end;


% LOOP OVER STATIONS, AVERAGING MAX AND MIN TO COMPUTE MEAN
for n = 1:m1;
   sd=Sd(n,:);
   sd=strtok(Sd(n,:),'.');
   sd=fliplr(deblank(fliplr(sd)));
   slet=sd(1);  % first letter
   snum=sd;
   snum(1)=[]; % number part of file in snum
   sd=['d' snum];
   se=['e' snum];
   
     
   eval(['Xd = load(''' pfd sd ''');']);
   D=Xd.X;
   yrD=D(:,1);
   D(:,1)=[];
      
   eval(['Xe = load(''' pfe se ''');']);
   E=Xe.X;
   yrE=E(:,1);
   E(:,1)=[];
   
   % Pull the common years of E and D
   yrgo = max([min(yrD) min(yrE)]);
   yrsp = min([max(yrD) max(yrE)]);
   Ld=yrD>=yrgo & yrD<=yrsp;
   Le=yrE>=yrgo & yrE<=yrsp;
   D=D(Ld,:);
   E=E(Le,:);
   yr1=(yrgo:yrsp)';
   
   % Average the max and min
   M=(D+E)/2;
   [mM,nM]=size(M);
   yrM=yr1;
   
   % Delete trailing and leading all-NaN years
   irow=(1:mM)';
   Lnone=  (all(isnan(M')))'; % cv pointing to all-NaN years
   if any(Lnone);
      irow(Lnone)=NaN;
      ivect=(nanmin(irow):nanmax(irow))';
      M=M(ivect,:);
      yrM=yr1(ivect);
   else;
   end;
   X=[yrM M];
   
   % Save file
   eval(['save ' pfm 'm' snum ' X;']);
   clear X M yrM;
   
end;
tspent=etime(clock,t0);
clc;
disp('v2conv3.m finished running');
disp(['   Number of files = ' int2str(m1)]);
disp(['   Time = ' num2str(tspent) ' seconds']);

   
            
      
        
   
  

   

   
   


   