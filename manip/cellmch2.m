function i=cellmch2(a,b)
% cellmch2:  
% L=cellmch1(a,b);
% Last revised 2-17-02
%
% Compares cell vectors a and b to find which members of a match members of b. Returns a column vector 
% i, same length as b, having either NaN or the index of the matched member in a.  Called by rwlset 
% in modifying the trimeach field of an rwlset
%
%*** INPUT
%
% a,b two cell vectors.  Elements are character strings.
%
%
%*** OUTPUT 
%
% i(?x1)i interger pointer to matched member of a; NaN if no match
%
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% Utility function for rwmeas/rwlspecs.  Looks at each member of b and finds the member (if any) of a that 
% matches it.  If a match, fills that slot of i with the index of the member in a.  If no match, inserts NaN.

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
i = repmat(NaN,nb,1); % initialize

for n = 1:nb; % loop over all series 
    i1 = strmatch(bb{n},aa,'exact');
    if isempty(i1);
        i(n)=NaN;
    elseif length(i1)>1;
        error([bb{n} ' matches two or more id names in the rwlset']);
    else;
        i(n)=i1;
    end;
end;



    
    
    
