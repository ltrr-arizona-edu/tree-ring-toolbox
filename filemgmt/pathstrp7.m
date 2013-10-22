function pathstrp7
% pathstrp7: strips user-specified paths from matlab path
% pathstrp7;
% Last revised 2004-10-22
%
% Strips Meko directories from MATLAB path.  
%
%*** INPUT --none
%*** OUTPUT-- none
%*** REFERENCES -- none
%*** UW FUNCTIONS CALLED -- none
%*** TOOLBOXES NEEDED -- none
%
%*** NOTES
%
% Used in testing GEOS 595E scripts to make sure I use only the versions in the directory sent to class
% If you edit this function, also edit pathfull.m
% 
% \mlb\filemgmt\ is not stripped because that is where pathstrp.m and pathfull.m are stored


dira='c:\mlb';
%dirb='c:\mlb\filemgmt';
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
rmpath(dira,dirc,dird,dire,dirf,dirg,dirh,diri,dirj,dirk,dirl,dirm,dirn,diro,dirp);

