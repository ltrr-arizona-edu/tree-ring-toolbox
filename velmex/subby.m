function rwlset=subby(rwlset)
% subby: utility subfunction to update obsolete structure of an rwlset with trimmall and trimeach fields
if isempty (rwlset.name{1}); % no rwlsets yet assigned
    rwlset.trimall{1}=[];
    rwlset.trimeach{1}=[];
else; % have at least one rwlset already; make its trimeach and trimall null
    nset = length(rwlset.name); % number of existing rwlsets
    for n = 1:nset; % loop over existing rwlsets;
        nser = length (rwlset.idnames{n});
        rwlset.trimall{n} = []; % 
        rwlset.trimeach{n}=cell(1,nser);
    end;
end;
