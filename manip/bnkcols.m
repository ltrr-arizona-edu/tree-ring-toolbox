function Bnk=bnkcols(nrows)
% bnkcols: utility function to build blank columns for spacing in tables
% Bnk=bnkcols(nrows);
% Last revised 9-24-02
%
% Utility function to build blank columns for spacing in tables.  Structure Bnk returned
% Bnk.n1 a blank char matrix of dimension nrows x 1,  Bnk.n2 a blank char matrix of 
% dimension nrows x 2, .... , out to dimension nrows x 9.  
%
%*** INPUT
%
% nrow (1 x 1)i   number of rows in desired char matrices
%
%*** OUTPUT
%
%  Bnk.n?  structure with the fields holding the blank char matrices
%   .n1  -- char mtx nrows x 1
%   .n2 --  blank char mtx nrows 02
%   :
%   .n9 -- blank ... nrows x 9
%
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES

% Hard code
maxblank=9;

% Make the structure of blank char matrices
for n=1:maxblank;
        eval(['Bnk.n' int2str(n) '= repmat(     blanks(' int2str(n) '),' int2str(nrows) ',1       );']);
end;