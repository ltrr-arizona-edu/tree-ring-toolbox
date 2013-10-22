function pathfull7
% pathfull7: builds full meko matlab path
% pathfull;
% Last revised 2006-12-14
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

dira='c:\mlb\filemgmt';
dirb='c:\mlb\velmex';
dirc='c:\mlb';
dird='c:\mlb\rwlook';
dire='c:\mlb\stz';
dirf='c:\mlb\convrt';
dirg='c:\mlb\manip';
dirh='c:\mlb\flowrec';
diri='c:\mlb\daily';
dirj='c:\mlb\watbal';
dirk='c:\mlb\screen';
dirl='c:\mlb\evol';
dirm='c:\mlb\sohel';
dirn='c:\mlb\mapdmm';
diro='c:\mlb\public';
dirp='c:\mlb\spatrec';


part1 = [' ' dira ' ' dirb ' ' dirc ' ' dird ' '];
part2 = [' ' dire ' ' dirf ' ' dirg ' ' dirh ' '];
part3 = [' ' diri ' ' dirj ' ' dirk ' ' dirl ' '];
part4 = [' ' dirm ' ' dirn ' ' diro ' ' dirp ' '];
pathy = [part1 part2 part3 part4];

eval(['addpath ' pathy ' -begin']);
