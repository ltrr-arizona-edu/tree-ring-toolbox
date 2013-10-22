function asort=alphtize(a)
% alphtize:  alphabetize cell strings
% CALL: aphatize;
%
% Meko 2-8-99
%
%******* INPUT
%
% a {1 x na}c  cell variable with na row strings
%
%******** OUTPUT
%
% asort{1 x na}c  cell variable with na row strings, in alphabetical order



if ~iscell(a);
   error('a must be cell');
end

[ma,na]=size(a);
if ma~=1;
   error('row-size of a must be 1');
end
% na is number of strings in cell variable a

% Storage for sorted cell
asort = cell(1,na);

% Convert a to a string matri
c = char(a);

% Sort c 
d = sortrows(c);

% Convert the char matrix to cells
for n = 1:na;
   atemp = strtok(d(n,:)); % lop off right blanks
   asort{n} = atemp; 
end

   