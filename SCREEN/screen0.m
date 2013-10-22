function screen0
%
% Make site-to-station distance matrix (km) and save it with
% related info for later use in screen1.m

% Get long-lat coords of tree sites
txt1='Input tree coords (map units) -- like c:\wrk6\xytreeA.mat';
txt2='Input clim coords (map units) -- like xyP.mat';
txt3='Input screen1.m misc setting -- like AP1IS1.mat';
[file1,path1] = uigetfile('xytree?.mat',txt1);
[file2,path2] = uigetfile('xyP.mat',txt2 );

char=input('Letter code for this tree-ring matrix: ','s');

fn3=[char '*.mat'];

[file3,path3] =  uigetfile(fn3,txt3);

eval(['load ',path1,file1])
eval(['load ',path2,file2])
eval(['load ',path3,file3])
k1=[1 2]; % hard code this dec2dms.m

% gcdist.m required deg,min,sec units
PT= dec2dms2(TC,k1); % change coordinates from map long-lat units
PC= dec2dms2(CC,k1); % likewise for clim station coords


disp('Sit back and relax; this could take 4 minutes')
P=gcdist(PT,PC); % calc dist (km) from each tree-site to each station
maxs=mx(1); % maximum number of stations considered
[D,W,N,I2]=near3(P,ds,maxs,maxs,1);
% D, W, N defined in gcdist.m;  I2 is pointer to valid chronologies, 
%		in other word, chronologies with at least one climate station
%		in the search radius
%
txt4='Output file name for distance info -- like dista.mat'
file4=uiputfile('DIST?*.mat',txt4);
eval(['save ',file4,' D W N I2'])


