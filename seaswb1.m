function seaswb1
% seaswb1:  seasonalize monthly water-balance quantities
% CALL: seaswb1
%
% Meko 10-7-97
%
%**********************  IN *****************************************
%
% No input arguments
%
% User prompted for name of file with following information on lines
%  1-number of sites
%  2-path to input .mat files with cell arrays  (holding the 16 variables)
%  3-begmo, endmo -- the start and end months of season 
%  4-yrgo, yrsp -- start and end year for the output matices
%  5-five letter season code: eg:  nv-ap
%
%
%*************************  OUT *************************************
%
% a 17-col tsm matrix is added to the site's input .mat file
%   col 1 -- year
%   col 2 -- pcp
%   col 3 -- tmp
%   col 4 -- RO  runoff
%   col 5 -- Z-index
%   col 6 -- X  PDSI
%   col 7 -- W  soil moisture (mid month)
%   cols 8-17 SI  stress index, which is P-cPE, where c = .1 ,.2, ... 1.0
%
% The name of the 17-col matrix is the 5-letter season code (see input)
%
%**********************************************************************

% Get path info and season specs
[file1,path1]=uigetfile('seas1.txt','Input control file');
pf1=[path1 file1];
fid1=fopen(pf1,'r');


% Number of sites
nsites = str2num(fgetl(fid1));

% Path to input files
path2=strtok(fgetl(fid1));
path2=fliplr(path2);
path2=fliplr(strtok(path2));
if ~strcmp('\',path2(length(path2)));
   error('Path does not end in \');
end



%----------------- Start, end month of season
ctemp = fgetl(fid1);
% get rid of leading and trailing blanks
ctemp=deblank(ctemp);
ctemp=deblank(fliplr(ctemp));
ctemp=fliplr(ctemp);
% take start as month ended by first blank
begmo = str2num(strtok(ctemp));
% cut off chars thru first space
i1 = find(isspace(ctemp));
ctemp(1:i1(1))=[];
% get rid of leading blanks
ctemp=fliplr(deblank(fliplr(ctemp)));
endmo = str2num(ctemp);



%----------------- Start, end year of master 17-col storage matrices
ctemp = fgetl(fid1);
% get rid of leading and trailing blanks
ctemp=deblank(ctemp);
ctemp=deblank(fliplr(ctemp));
ctemp=fliplr(ctemp);
% take start as month ended by first blank
yrgo = str2num(strtok(ctemp));
% cut off chars thru first space
i1 = find(isspace(ctemp));
ctemp(1:i1(1))=[];
% get rid of leading blanks
ctemp=fliplr(deblank(fliplr(ctemp)));
yrsp = str2num(ctemp);


%--------------  5-letter season code
seascode = fgetl(fid1);
seascode = strtok(seascode);



%*****************  LOOP OVER STATIONS   ***************************
a=NaN;
nyrs = yrsp-yrgo+1;

for n = 1: nsites;
   
   % Allocate for the 17-col matrix
   XX = a(ones(nyrs,1),ones(17,1));
   yrbig = (yrgo:yrsp)';
   XX(:,1)=yrbig;
   
   
   % Load file with the pdsi1.m output data
   pfin = [path2 'wbout' num2str(n)];
   eval(['load ' pfin]);
   
   %------------ Seasonalize various data
   
   % pcp
   
   A=datout{7};
   % find first, last year of monthly data
   yr3 = A(:,1);
   yrs = [min(yr3) max(yr3)];
   F = seaspt(A,begmo,endmo,yrs,1);
   % find first, last year of seasonalized data
   yr4=F(:,1);
   yr4go=min(yr4);
   yr4sp=max(yr4);
   % Pointer to storage matrix
   L1 = yrbig>=yr4go & yrbig<=yr4sp; 
   % store seasonalized data
   XX(L1,2)=F(:,2);
   
 % tmp 
   A=datout{8};
   % find first, last year of monthly data
   yr3 = A(:,1);
   yrs = [min(yr3) max(yr3)];
   F = seaspt(A,begmo,endmo,yrs,2);
   % find first, last year of seasonalized data
   yr4=F(:,1);
   yr4go=min(yr4);
   yr4sp=max(yr4);
   % Pointer to storage matrix
   L1 = yrbig>=yr4go & yrbig<=yr4sp; 
   % store seasonalized data
   XX(L1,3)=F(:,2);
 
   % model runoff 
   A=datout{5};
   % find first, last year of monthly data
   yr3 = A(:,1);
   yrs = [min(yr3) max(yr3)];
   F = seaspt(A,begmo,endmo,yrs,1);
   % find first, last year of seasonalized data
   yr4=F(:,1);
   yr4go=min(yr4);
   yr4sp=max(yr4);
   % Pointer to storage matrix
   L1 = yrbig>=yr4go & yrbig<=yr4sp; 
   % store seasonalized data
   XX(L1,4)=F(:,2);
   
     % Z-index 
   A=datout{1};
   % find first, last year of monthly data
   yr3 = A(:,1);
   yrs = [min(yr3) max(yr3)];
   F = seaspt(A,begmo,endmo,yrs,2);
   % find first, last year of seasonalized data
   yr4=F(:,1);
   yr4go=min(yr4);
   yr4sp=max(yr4);
   % Pointer to storage matrix
   L1 = yrbig>=yr4go & yrbig<=yr4sp; 
   % store seasonalized data
   XX(L1,5)=F(:,2);
 
   % X = pdsi 
   A=datout{2};
   % find first, last year of monthly data
   yr3 = A(:,1);
   yrs = [min(yr3) max(yr3)];
   F = seaspt(A,begmo,endmo,yrs,2);
   % find first, last year of seasonalized data
   yr4=F(:,1);
   yr4go=min(yr4);
   yr4sp=max(yr4);
   % Pointer to storage matrix
   L1 = yrbig>=yr4go & yrbig<=yr4sp; 
   % store seasonalized data
   XX(L1,6)=F(:,2);
   
   % W = ave soil moisture (mid-month) 
   A=datout{4};
   % find first, last year of monthly data
   yr3 = A(:,1);
   yrs = [min(yr3) max(yr3)];
   F = seaspt(A,begmo,endmo,yrs,2);
   % find first, last year of seasonalized data
   yr4=F(:,1);
   yr4go=min(yr4);
   yr4sp=max(yr4);
   % Pointer to storage matrix
   L1 = yrbig>=yr4go & yrbig<=yr4sp; 
   % store seasonalized data
   XX(L1,7)=F(:,2);
   
   % Now the adjusted pcp:  pcp minus a fraction of PE, where fraction is
   % .1, .2, ... 1.0
   A=datout{6};
   % find first, last year of monthly data
   yr3 = A(:,1);
   yrs = [min(yr3) max(yr3)];
   
   % Loop over the fractions
   for j = 1:10;
      B = A(:,:,j);
      F = seaspt(B,begmo,endmo,yrs,1);
      % find first, last year of seasonalized data
      yr4=F(:,1);
      yr4go=min(yr4);
      yr4sp=max(yr4);
      % Pointer to storage matrix
      L1 = yrbig>=yr4go & yrbig<=yr4sp; 
      % store seasonalized data
      XX(L1,7+j)=F(:,2);
   end
   
   % Rename XX according to the season code
   eval([seascode ' = XX;']);
   
   
   % Add XX to the existing watbalm1.m output file
   eval(['save ' pfin ' ' seascode ' -append;']); 
   
    
end
