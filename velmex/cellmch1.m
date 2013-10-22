function L=cellmch1(a,b)
% cellmch1:  compares cell vectors a and b and returns 1 for cells of b matching any cell in a
% L=cellmch1(a,b);
% Last revised 2-17-02
%
% compares cell vectors a and b and returns 1 for cells of b matching any cell in a. Called by
% rwlspecs.m in modifying an rwlset
%
%*** INPUT
%
% a,c two cell vectors.  Elements are character strings.
%
%
%*** OUTPUT 
%
% L(?x1)L logical variable marking cells in b that match any in a
%
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% Utility function for rwmeas/rwlspecs

if ~(iscell(a) & iscell(b));
    error('a and b must be cells');
end;
[ma,na]=size(a);
[mb,nb]=size(b);
if ~(ma==1 & mb==1);
    error('Row dim of a and b must be 1');
end;

aa=lower(a'); % make into col-cell, and lower case for comparison
bb=lower(b'); % make into col-cell
L = logical(zeros(nb,1)); % initialize

for n = 1:na; % loop over all series currently in the rwlset
    i1 = strmatch(aa{n},bb,'exact');
    L(i1)=1;
end;



    
    
    
