function pathfull
% pathfull: builds full meko matlab path
% pathfull;
% Last revised 8-23-01
%
% Adds full set of Meko directories to MATLAB path.  
%
%*** INPUT --none
%*** OUTPUT-- none
%*** REFERENCES -- none
%*** UW FUNCTIONS CALLED -- none
%*** TOOLBOXES NEEDED -- none
%
%*** NOTES
%
% Typically, run this after upgrade to new version of MATLAB.  You will need to change cwd to
% c:\mlb\filemgmt\ beforehand so that pathfull.m is found. 
%
% If you edit this function, also edit pathstrp.m


dira='c:\mlb';
dirb='c:\mlb\filemgmt';
dirc='c:\mlb\rwlook';
dird='c:\mlb\stz';
dire='c:\mlb\convrt';
dirf='c:\mlb\manip';
dirg='c:\mlb\flowrec';
dirh='c:\mlb\daily';
diri='c:\mlb\watbal';
dirj='c:\mlb\screen';
dirk='c:\mlb\evol';
dirl='c:\mlb\sohel';
dirm='c:\mlb\mapdmm';
dirn='c:\mlb\public';
diro='c:\mlb\spatrec';
dirp='c:\mlb\velmex';

part1 = [' ' dira ' ' dirb ' ' dirc ' ' dird ' '];
part2 = [' ' dire ' ' dirf ' ' dirg ' ' dirh ' '];
part3 = [' ' diri ' ' dirj ' ' dirk ' ' dirl ' '];
part4 = [' ' dirm ' ' dirn ' ' diro ' ' dirp ' '];
pathy = [part1 part2 part3 part4];

eval(['addpath ' pathy ' -begin']);
