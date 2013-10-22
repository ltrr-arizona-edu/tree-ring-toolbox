function screen5(R1,R2,K2,iw,d,ny)
%
% make  ascii file of model statistics that can be
% used in surfer maps
%
%
% D Meko 6-28-96
%
% The 20 columns of the file to be produced are as follows:
%
% A x coordinate in dec deg longitude, negative west
% B y coordinate in dec deg lat
% C sequence number of tree-ring chronology, corresponding to
%   the column in X in d:\wrk6\tsm?.mat
%
%*** next columns deal with the model fit to full period using
%	  screen3.m
%
% D class of the model (1 or 2) based on site-station distance
%	 and completeness of data coverage
% E squared correlation coefficient u(t) and y(t);  equivalent
%	 to decimal proportion of tree-ring variance explainable
%	 with simple no-lag model of year-t ring on year-t climate
%	 Got by squaring R1 as ouput from screen3.m
% F decimal proportion of variance explained by the best-fit
%	 OE model.  S as output by screen3.m
% G order of B operator in the OE model
% H order of F operator in the OE model
% HH any significant OE model params at 2 std devs 1-yes, 0-no
% I number of significant (at 99%) autocorrelation coefs of
%	 residuals in lags 1-10 years
% J number of significant (at 99%) cross-correlations between
%   residuals and input in OE model
% K which climate station in the model; entry corresponds to 
%	 a column number of the climate matrix used by screen3.m
% L distance (km) from the tree-site to the climate station
% M number of years in modeling period
%
%*** next columns deal with the 5 models fit to full period
%  using screen1.m 
%
% N maximum proportion variance explained of the 5 station models
% O median  "
% P minimum "
% Q distance (km) from site to nearest of "5 nearest stations"
% R distance (km) from site to farthest of "5 nearest stations"
% S binary (1 or 0) flag telling that simple r (clim vs tree) 
%   is negative (1)
%**********************************************************

% Help in file selection by specifying letter code for tree-ring set
char=input('Letter code for tree-ring matrix: ','s');

% Get the plotting coordinates
fn5=['crnxy' char '.dat'];
txt5='Input file with tree site lats, lons';
[file5,path5]=uigetfile(fn5,txt5);
pf5=lower([path5 file5]);
eval(['load ' pf5]);  
set1=' crnxy';
fxy = ['crnxy' char];
eval(['AZ=' fxy '(:,1);']);
eval(['BZ=' fxy '(:,2);']);

% Get site-to-station distance info for nearest 5 stations
fn5=['dist' char '.mat'];
txt5='dist?.mat (in D:\wrk0)-- distance info';
[file5,path5]=uigetfile(fn5,txt5);
pf5=lower([path5 file5]);
eval(['load ' pf5]);  
QZ=(min(D'))';  % distance to closest of 5 nearest
RZ=(max(D'))';  % distance to fartherst of 5 nearest


% Get screen2.m output
fn5=[char '*.mat'];
txt5='?P#OS2.MAT (in D:\wrk0)-- screen2 output';
[file5,path5]=uigetfile(fn5,txt5);
pf5=lower([path5 file5]);
eval(['load ' pf5]);  

nsites=length(iw); % number of tree sites
CZ = (1:nsites)';
DZ=iw;
KZ =w;  % climate station number
LZ =d; % km station to site

% Get the screen3.m output
fn5=[char '*.mat'];
txt5='?P#OS3.MAT (in D:\wrk0)-- screen3 output';
[file5,path5]=uigetfile(fn5,txt5);
pf5=lower([path5 file5]);
eval(['load ' pf5]);  

EZ=R1(:,1) .^2; % squared r between u(t) and y(t)
SZ=R1(:,1)<0;  % flag for negative correlation, u(t) vs y(t)
FZ=S(:,1); % 1-(var(errors)/variance(y))
GZ=OB(:,1); % B-order
HZ=OF(:,1); % F-order
HHZ=K2(:,1); % 1=any model parms signif at 2 sd, 0=none
IZ=H(:,1); % number of cross-corrs (u vs e) out of 1st 10 sign at 99%
JZ=G(:,1); % number of autocorrs of errors out of 1st 10 ...
MZ=NY(:,1); % number of years final model based on

clear R1 R2 S OB OF H G NY

% Get the screen1.m output
fn5=[char '*.mat'];
txt5='?P#OS1.MAT (in D:\wrk0)-- screen1 output';
[file5,path5]=uigetfile(fn5,txt5);
pf5=lower([path5 file5]);
eval(['load ' pf5]);  

NZ=(max(S'))'; % of 5 nearest stations, highest variance explained fraction
OZ=(median(S'))'; % of ..., the median
PZ=(min(S'))';  % ... minimum...


ZZ=[AZ BZ CZ DZ EZ FZ GZ HZ HHZ IZ JZ KZ LZ MZ NZ OZ PZ QZ RZ SZ];

fmt1a='%7.2f %7.2f %4d %3d %5.2f %5.2f ';
fmt1b='%2d %2d %2d    %3d %3d   %4d %6.1f %4d ';
fmt1c='%5.2f %5.2f %5.2f %4.0f %4.0f %1d\n';
fmt1=[fmt1a fmt1b fmt1c];

% Open file for output ascii table
fn5=[char '*.dat'];
txt5='?P#OS5.DAT (in D:\wrk0)-- screeen5.m output ascii file';
[file5,path5]=uiputfile(fn5,txt5);
pf5=lower([path5 file5]);
ff=lower([path5 file5]);

fid1=fopen(ff,'w');

for n=1:nsites;
	zz=ZZ(n,:);
	fprintf(fid1,fmt1,zz);
end


fclose all
