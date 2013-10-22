function seqnum1
% seqnum1: add column of sequence numbers to 2-col decimal latitude - longitude file
% CALL: seqnum1
% USE: From foxpro and dms2dec.m, made a 2 col long-lat file for tree ring sites.
%   Need to plot sequential site number on surfer map, so need to slap on a column
%   of sequential numbers.  
%
% Meko 3-22-97
%
%******************* IN  (no args) *********************************************
%
% Files prompted for:
%   - .dat file of long and lat, in decimal degree units
%       Assumed to be 2-col numeric, and loadable into matlab with load command
%       
%
%**************** OUT (no args) ****************************************
%
% Files prompted for: 
%
%  - .dat filename to store the output 3-col matrix.  sequence number as col 3.
%      Three digits to right of decimal point in long, lat
%
%


% Get input file
[file1,path1]=uigetfile('*.dat','Input .dat file with lon, lats in dec degs');
pf1=[path1 file1];
eval(['load ' pf1]);

% Put raw long lat data in X
name1=lower(strtok(file1,'.'));
X=eval(name1);

% Size X and build augmented matrix
[mX,nX]=size(X);

j=(1:mX)';
Y=[X j];

% Format for output line
fmt='%8.3f %8.3f %4.0f\n';

% Open output file
[file2,path2]=uiputfile('*.dat','Output long-lat file with seq nos as col 3');
pf2=[path2 file2];
[fid2,message]=fopen(pf2,'w');
if fid2<0;
   disp(message)
else
   fprintf(fid2,fmt,Y');
end
fclose(fid2);




  



