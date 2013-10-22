function dispmon1
% Screen display of monthly climate data. Useful for quality control in updating.
% dispmon1;
% Last revised 1-21-00
%
% Write monthly precip data for selected years on the screen. 
%
%
%*** INPUT
%
% No args.
% User prompted for information.
%
%
%*** OUTPUT 
%
% None
%
%*** REFERENCES
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
% NOTES
%
% Review "hard code" settings.  Set these to best automate work.

% Hard Coded 
fmt1='%4d %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n';
fscale=100 % ppt from inches to hundredths
X='Z'; % matrix holding monthly climate data in input
yrgo=1994; % first year of display period
yrsp=1999; % last year ...

% Get file
[file1,path1]=uigetfile('*.mat','Input .mat file of climate data');
pf1=[path1 file1];
eval(['load ' pf1 ';']);
eval(['X=' X ';']);

% Pull desired years
L = X(:,1)>=yrgo & X(:,1)<=yrsp;
if ~any(L);
   error(['No years between ' int2str(yrgo) ' and ' int2str(yrsp)]);
else;
   X=X(L,:);
end;

% Scale
X(:,2:13)=X(:,2:13)*fscale;
sprintf(fmt1,X')




