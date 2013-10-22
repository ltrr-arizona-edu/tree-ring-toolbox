function eflag=corenm1(id)
% corenm1: check syntax of core ids
% corenm1(id);
% Last revised 10-10-01
%
% Checks that core ids in cell array all follow naming convention that:
%   Begin with 3-char site code, all letters
%   Next comes numeric tree number
%   Next comes core letter
%   Optionally, next comes sequence number within core
%
%*** INPUT
%
% id{} cell array of core ids <{'pdf02a','pdf03b'}> (See notes)
%
%*** OUTPUT 
%
% No arguments.
% Error message if any invalid ids
%
%*** REFERENCES --- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES


nid=length(id);

for n =1:nid;
    nameser=id{n};
    
    % Length of name
    nlen = length(nameser);
    if nlen>8;
        error(['Series ID ' nameser ' has more than 8 chars']);
    elseif nlen<5;;
        error(['Series ID ' nameser ' needs at least 5 chars for a site code, tree no and core letter']);
        
    end;
    
    
    % First 3 chars must be letters
    d1 = nameser(1:3);
    if ~all(isletter(d1));
        error([nameser ' first 3 chars must be letters']);
    end;
    sitecode=d1;
    clear d1;
    
    % Check that have a core letter, following a tree number
    d1=nameser;
    d1(1:3) = []; % strip off site code
    i1 = find(isletter(d1)); % find remaining letters
    if isempty(i1);
        error([nameser ' has no letters after the initial -- needs a core letter']);
    else;
        treeno=d1(1:i1(1)-1);
        if length(treeno)==1;
            treeno=['0' treeno];
        elseif length(treeno)>3;
            error([nameser ' tree number should be max of 999 and max of 3 digits']);
        end;
    end;
    
    % Check core letter and for trailing core-segment number
    d1 = d1(i1(1):length(d1));
    i1=find(isletter(d1));
    if length(d1)==1;
        coreletter=d1;
        nameser=[sitecode treeno coreletter];
    else;
        % Strip the letter off, and should have only a 1-digit number
        coreletter=d1(i1(1));
        d1(1)=[];
        if length(d1)>1;
            error([nameser ' has a core-segment number after the core letter, but this should be length 1']);
        else;
            if isletter(d1); 
                error(['Last char of ' nameser ' should be a number for this case']);
            end;
            segno = d1;
            nameser=[sitecode treeno coreletter segno];
        end;
    end;
    
end; % end of loop over ids
disp('corenm1 ran OK');
