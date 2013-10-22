function y=uniqcell(a,b)
% uniqcell:  compares two cell arrays of names and returns unique (non-duplicate) names
% y=uniqcell(a,b);
% Last revised 9-27-01
%
% Find the unique names in two cell arrays of names
%
%*** INPUT
%
% a,c two string arrays.  Elements are character strings.
%
%
%*** OUTPUT 
%
% y{} cell array of names from a and b, with duplicates omitted
%
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% Written as utility for rwc2rwl.m to identify core ids from lists of core ids in total ringwidth and partial
% ringwidth files
%
if (isempty(a{1}) & isempty(b{1}));
    y=[];
elseif isempty(a{1});
    y=b;
elseif isempty(b{1});
    y=a;
else;
    c=[a b];
    nc=size(c,2); % number of series
    d=c;
    
    % compare elements of c and d
    c0=c;
    j=0;
    for n =1:nc;
        cthis = c0{n};
        L=strcmpi(cthis,d);
        if any(L);
            d(L)=[];
            j=j+1;
            g{j}=cthis;
            if isempty(d);
                y=g;
                return;
                
            end;
        end;
    end;
end;

        
        
        
        
    
    
    
