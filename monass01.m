function [nsites,L,t,y]=monass01
% monass01:  assess completeness of monthly data file for specified subperiod
%
%
%**************** IN ****************
%
% No input args
%
% User prompted to select input file with following info
%
% line 1 -- path to input .mat files (e.g., c:\data\infiles\)
% line 2 -- Length of reference period
% line 3 --  first end yr, last end yr 
% lines 4 ...  prefixes of .mat files holding monthly data


%*************** OUT ***********
%
% nfiles (1 x 1)i number of input files
% L (? x 1)L  logical vector telling which of the series passes sufficiency test

%-------------- GET INPUT CONTROL AND READ FIRST 3 LINES
close all

[file1,path1]=uigetfile('*.txt','input control file');
pf1=[path1 file1];
fid1=fopen(pf1,'r');

% Read path to input .mat files
path2=fgetl(fid1);
path2=strtok(path2);
path2=strtok(fliplr(path2));
path2=fliplr(path2);

% Read length of analysis period
c=fgetl(fid1);
nlength = str2num(strtok(c));


% Get analysis period end yrs
c=fgetl(fid1);

c1=strtok(c);
c1=strtok(fliplr(c1));
c1=fliplr(c1);
yrgo=str2num(c1);

c2=strtok(fliplr(c));
c2=strtok(fliplr(c2));
yrsp=str2num(c2);

% Get required minimum pctg of non-NaN data for any column
c = fgetl(fid1);
pcrit = str2num(strtok(c));


% Compute number of analysis periods
nperiods = yrsp-yrgo+1;


%********* LOOP OVER SITES

clc
disp(['Counting Number of Files']);
nsites=0;
kwhile=1;
while kwhile==1;
   c=fgetl(fid1);
   if ~(feof(fid1) | length(c)<2);
      nsites=nsites+1;
   else
      kwhile=0;
   end
end
disp(['   Done: Number of files is ' int2str(nsites)]);
disp(['         Analysis period: ' int2str(yrgo) '-' int2str(yrsp)]);


%************* Prompt for name of matrix holding data in each .mat file

 prompt={'Enter name of data matrix:'};
 def={'X'};
 tit='Name of Matrix Holding Data in .mat Files';
 lineNo=1;
 varnm=inputdlg(prompt,tit,lineNo,def);
 varnm=varnm{1};

%***********  LOOP OVER SITES AGAIN

a = NaN;

% Size logical matrix to hold indicator whether station passes or not
L =repmat(a,nsites,nperiods);

frewind(fid1);
c=fgetl(fid1);
c=fgetl(fid1);
c=fgetl(fid1);
c=fgetl(fid1);
nyrs = yrsp-yrgo+1;

for n =1:nsites
   flpref = strtok(fgetl(fid1)); % prefix of .mat file
   pf2=[path2 flpref];
   eval(['load ' pf2]);
   eval(['X = ' varnm ';']);
   yr = X(:,1);
   X=X(:,2:13);
   
   for m = 1:nperiods
      yr1 = yrgo+m-1;
      yr2 = yr1 + nlength -1;
      L1 = yr>=yr1 & yr<=yr2;
      if any(L1)
         X1=X(L1,:);
         sum1 = sum(~isnan(L1));
         if (sum1/nlength) < pcrit;
            L(n,m)=0;
         else
            L(n,m)=1;
         end
      else
         L(n,m)=0;
      end
      
   end
end

L=logical(L);
  
fclose all

% Plot result
figure(1);
t = (yrgo:yrsp)'; % end year of period
y = (sum(L))'; % number of stations passing test
plot(t,y);
grid;
zoom xon
txt1=['Number of stations passing out of ' int2str(nsites)];
title(txt1);
ylabel('Number of stations');
txt2=['Ending year of ' int2str(nlength) '-yr period']; 
xlabel(txt2);