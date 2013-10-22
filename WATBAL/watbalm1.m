function watbalm1
% watbalm1:  compute water-balance quantities for multiple site-centered climate
% CALL: watbalm1
%
%*********************  IN ***********************************
%
% No input arguments
%
% A control file prompted for, with following info line by line
%
%   The number of sites <57>
%   start, end years of desired watbal output, assumed same for each site
%   start, end years of 'CAFEC' period, "
%   prefix to the site-centered series <grdm>, assumed same for pcp and tmp
%   path\file to file with names of tree sites <c:\...\wrktree\crnfna.txt>
%   path\file to long-lat file for sites <c:\...\wrktree\crnxya.dat>
%   path\file to one-line or many-line soil moisture specs 
%      <c:\...\wrktree\soilmois.dat>
%   path to the site-centered input monthly pcp files
%      <c:\...\pcpfiles\regout\>
%   path to the site-centered input monthly tmp files
%      <c:\...\tmpfiles\regout\>
%   path\file to lookup table dayz <c:\mlb\watbal\daylennh.mat>
%   path\file to .mat file with cell arrays snowinf, kopt, penopts, datpen
%      <c:\...\wrktree\pdmisc.mat>
%   path to directory to hold the watbalm1.m output
%      <c:\projs\af1\sacrflow\outwat\>  (ideally empty)
%
%
%*****************   OUT ************************************
%
% In the output directory specified above, a .mat file for each
% site, holding the site's watbalm1.m output.  'wb' starts the file name,
% and the sequential site number ends it. Example:  wb1.mat, wb2.mat,...
% wb57.mat.  
%
% Each .mat file holds a cell array datout with the 6 variables:
% Z,X,XM,W,RO,SI, as defined in pdsi1.m
%
%


% Get master control file
[file1,path1]=uigetfile('*.txt','Control file');
pf1=[path1 file1];
fid1 = fopen(pf1,'r');

nsites=str2num(fgetl(fid1)); % Number of sites


%---------------- Get start end year of desired output period
ctemp = fgetl(fid1);
% get rid of leading and trailing blanks
ctemp=deblank(ctemp);
ctemp=deblank(fliplr(ctemp));
ctemp=fliplr(ctemp);
% take start as year ended by first blank
yr1go = str2num(strtok(ctemp));
% cut off chars thru first space
i1 = find(isspace(ctemp));
ctemp(1:i1(1))=[];
% get rid of leading blanks
ctemp=fliplr(deblank(fliplr(ctemp)));
yr1sp = str2num(ctemp);

%---------------- Get start end year of cafec period
ctemp = fgetl(fid1);
% get rid of leading and trailing blanks
ctemp=deblank(ctemp);
ctemp=deblank(fliplr(ctemp));
ctemp=fliplr(ctemp);
% take start as year ended by first blank
yr2go = str2num(strtok(ctemp));
% cut off chars thru first space
i1 = find(isspace(ctemp));
ctemp(1:i1(1))=[];
% get rid of leading blanks
ctemp=fliplr(deblank(fliplr(ctemp)));
yr2sp = str2num(ctemp);

%-------------  Get prefix to site-centered pcp, tmp files
prefix1=strtok(fgetl(fid1));



%----------------- Matrix of names of tree sites
blnks20 = blanks(20);
sitenms = repmat(blnks20,nsites,1);
pf2= fgetl(fid1);  % name of .dat file holding names of .crn files

% Open the file of filenames 
fid2=fopen(pf2,'r');

% Build the matrix of names
for n = 1:nsites;
   ctemp = fgetl(fid2);
   ctemp = fliplr(deblank(fliplr(ctemp)));
   ctemp = strtok(ctemp,'.');
   nlen1 = length(ctemp);
   sitenms(n,1:nlen1)=ctemp;
   temp1 = [prefix1 num2str(n)];
   nlen2=length(temp1);
   sitenms(n,nlen1+1)='-';
   igo1=nlen1+2;
   isp1 = igo1 + nlen2-1;
   sitenms(n,igo1:isp1)=temp1;
end
fclose(fid2);


%-------------- Longitudes and latitudes
pf3=fgetl(fid1);
eval(['load ' pf3 ]); % pf3 is path\file for file of long and lat
file3=strtok(fliplr(strtok(fliplr(pf3),'\')),'.'); % prefix of lon-lat filename file
eval(['lonlat = ' file3]);


%------------- Soil moisture in top and underlying level, inches
%
% File with data can be 1 row, meaning all sites use same values, or can be nsites rows
pf4 = fgetl(fid1); 
eval(['load ' pf4 ]); % pf4 is path\file for soil moisture
file4=strtok(fliplr(strtok(fliplr(pf4),'\')),'.'); % prefix of soilmoisture file
eval(['awc = ' file4]);
if size(awc,1)==1; % all sites use same awc values
   ksoil=1; % use single value of awc for all sites
else;
   ksoil=2; % vary awc by site
end


%--------------  Get Paths for monthly pcp and tmp files
pathpcp=fgetl(fid1);
pathtmp=fgetl(fid1);


%----------------- Get daylength factor table
pf5=fgetl(fid1);
eval(['load ' pf5]);  % dayz should be table

%---------------- Get miscellaneous pdsi control
pf6=fgetl(fid1);
eval(['load ' pf6]);  % gives snowinf,datpen,penopts,kopt

%------------------- Get dir to put output in
path7 = strtok(fgetl(fid1));
len1 = length(path7);
if ~strcmp(path7(len1),'\');
   error('Directory for output must end with \');
end


%****************  SITE-BY-SITE MODELING *****************************


 % First set variables that remain same from site to site
 
 % Available water capacity
 if ksoil==1; % if using same awc settings for all sites
    watcap=[awc(1) awc(2)];
 end
 
 % Years for output and for cafec period
 YRS = [yr1go yr1sp; yr2go yr2sp];
 
 
 for n = 1:nsites; % loop over sites
    
    % Build output file name
    pf7= [path7 'wbout' int2str(n)];
    
    
    % Get file prefix for pcp and tmp file
    pf8=[pathpcp prefix1 int2str(n)];
    pf9=[pathtmp prefix1 int2str(n)];
    
    % build datmon
    clear Z; % just in case
    eval(['load ' pf8]); % pcp data in Z
    P = Z; 
    clear Z;
    eval(['load ' pf9]); % tmp data in Z
    T=Z;
    clear Z;
    datmon = {P,T};
    
    % build datother
    lat = lonlat(n,2); % latitude, decimal degrees
    if ksoil==2;
       watcap=[awc(n,1) awc(n,2)];
    end    
    datother={lat,watcap,YRS,dayz};
    
    
    % Text input for plot titles
    textin = {sitenms(n,1:20)};
    
    
    % Run the water balance    
    datout=pdsi1(datmon,datother,snowinf,textin,kopt,penopts,datpen);
    
    
    % Store output 
    eval(['save ' pf7  ' datout']);
 end
 

fclose all;  


