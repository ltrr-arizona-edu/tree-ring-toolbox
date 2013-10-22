function b=deblankc(a);
% deblankc:  deblank the elements (rows) of a col-cell of strings
% b=deblankc(a);
% Last revised 5-4-99
%
% Used in preparing Santa Rita precip station labels.  After using edix
% and matlab to make a col-cell of strings for station labels, saw that
% each cell was padded to length 5 with blanks.  Needed to remove the
% blanks
%
%*** IN 
%
% a {} col-cell of strings, some right-padded with blanks
%
%
%*** OUT 
%
% b {} col-cell of strings, deblanked
%
%*** REFERENCES -- none
%*** UW FUNCTIONS CALLED -- none
%
%*** NOTES

b=a; % allocate for revised cell
[m1,n1]=size(a);

% loop over cell elements
for n = 1:m1;
   b{n}=deblank(a{n});
end



