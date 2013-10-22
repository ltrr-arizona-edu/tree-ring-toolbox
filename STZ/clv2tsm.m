function tsm = clv2tsm
% strung-out column vector to time series matrix

clear all

% Prompt for the .mat file name

% Get the .mat filename interactively; load file
flmat=uigetfile('*.mat','.MAT filename ?');
flold=flmat; % Will need this file name later
eval(['load ',flmat]);
clear flmat;


xtemp = input('Strung-out vector name: ','s');
nametemp=input('Names Label string matrix: ','s');
yrstemp=input('Years and row index matrix: ','s');

outfname = input('Desired name of output tsm: ','s');


% Find Inclusive period
eval(['years= ',yrstemp,';']); 
eval(['names= ',nametemp,';']);% 
eval(['xcv = ',xtemp,';']); % strung-out col vector
yfirst = min(years(:,1));
ylast = max(years(:,2));


yrvect= (min(years(:,1)):max(years(:,2)))';

% Compute target row indices in times series matrix
rowtsm = [(years(:,1)-yfirst+1) (years(:,2)-yfirst+1)];

% Get the individual series
[ns,dum1]=size(years);
mout = ylast-yfirst+1;
nout = ns+1;
a=NaN;
tsm = a(ones(mout,1),ones(nout,1));
tsm(:,1)=(yfirst:ylast)';

rowin  = [years(:,3)  years(:,3)+years(:,2)-years(:,1)];
	
for k = 1:ns
	rows1 = rowin(k,:);
	rows2 = rowtsm(k,:);
	tsm(rows2(1):rows2(2),k+1) = xcv(rows1(1):rows1(2));
end
