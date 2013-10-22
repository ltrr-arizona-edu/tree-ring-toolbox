% ascrec1.m  ascii file of reconstructed and actual precip from
% spatrec1.m output.  For Chris Baisan, 1-8-97

% Change to the data directory with the .mat file
cd c:\wrk9

% Load the desired reconstruction and related info
load that8e;  
% reconstruction vector is in YP, Observed calib ppt is in Y;
% corresp year vectors are yr4 and yr2

% Plot to check out the time series
plot(yr2,Y,yr4,YP)


% Open an output file for ascii output
[file1,path1]=uiputfile('*.dat','Ascii output file');
pf1=[path1 file1];
fid1=fopen(pf1,'w');

% Set up format for output
fmt1='%4.0f %8.2f\n';

% Splice the reconstruction (1271-1915) to the observed (1916-68)
z1=[yr2 Y];
z2=[yr4 YP];
z3=[z2;z1];


% Write it
fprintf(fid1,fmt1,z3');

fclose all


