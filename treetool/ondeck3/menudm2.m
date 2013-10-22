function L=menudm2(tit,names)
% menudm2: make a menu and logical pointer to selected item
% L=menudm2(tit,names);
% Written 5-19-98; Last revised 8-26-99
%
%*****************  IN
%
% names {}s    1x? cell matrix of names in menu
% tit (1 x ?)s header for menu items (for call to menu.m)
%
%****************** OUT
%
% L (1 x ?)L  pointer to selected menu items: selected(1) or not selected (0)
%
%*** NOTES
%
% Subfunction of crvfit.m that changes '-f'  (already fit) to '-n' (not fit)
% ALSO USED BY TREEI


ntot = size(names,2); % number of items to pick from
nchar = size(names{1},2); % number of characters in each name, including pre and suffix

% Initialize logical pointer to say none yet picked
L=logical(zeros(ntot,1));

% Make screen menu
kmen1 = menu(tit,names);

% Change the last 2 characters of the name of the selected item
names{kmen1}((nchar-1):nchar)='-n';
L(kmen1)=1;


      