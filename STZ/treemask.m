function tmask=treemask(Tnms,tmask)
% treemask:  build or revise mask for omitting trees from chronology 
% CALL: tmask=treemask(Tnms,tmask);
%**************
%
%************** IN
%
% nms (? x 8)s names of cores or trees
% tmask (? x 1)L  mask (0==mark for deletion, 0==keep)
%
%***************** OUT
%
% tmask (? x 1)L mask: 
%   ==0 mark tree for omission
%   ==1 keep tree
%
%************* NOTES
%
% User decides on masking from statistics and knowledge of data

% Compute number of trees.  This would have been set by treei.m, which considers 
% the core mask
[ns,ms]=size(Tnms); % ns is number of trees

% Fill the string matrix indicating whether tree already marked for excluding
strlive='-n';
strdie='-m';

% Initialize string matrix telling whether tree to be kept or excluded
FV=repmat(strlive,ns,1);
% If cmask does not already exist, set to not mask any cores
if isempty(tmask) | ~(exist('tmask')==1);
   tmask=logical(ones(ns,1));
else; % if tmask does exist, set the mask
   if any(~tmask); % if any trees previously marked for excluding
      nexclude = sum(~tmask);
      FV(~tmask,:)=repmat(strdie,nexclude,1);
   end;
end;



% --- Build the mask

% Build string vector for "Quit" in a menu
dum=[];
for i=1:ms-4,
  dum=[dum,' '];
end
qt=['QUIT',dum];


k=0;
tmask=ones(ns,1);
while k<ns+1,
  nmsn=[Tnms FV;qt '  '];
  k=find(menudm2('Change tree mask?',(cellstr(nmsn))'));
  if k<=ns; % change the tree mask for this tree
     if tmask(k);
        tmask(k)=0;
        FV(k,:)=strdie;
     else;
        tmask(k)=1;
        FV(k,:)=strlive;
     end;
  end;
end;

  

% Make mask logical
tmask=logical(tmask);
