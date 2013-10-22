function cmask=coremask(nms,cmask)
% coremask: builds or revises the mask for excluding cores from chronology
% CALL: cmask=coremask(nms,cmask)
%
%************** IN
%
% nms (? x 8)s names of cores or trees
% cmask (? x 1)L  mask (0==exclude, 1==include)
%
%*** UW FUNCTIONS CALLED
%
% menudm2
%
%************* NOTES 
%
% User bases decision to mask on statistics and other information. 
% Only unmasked series (cmask==1)will be used in subsequent analyses.

% Find out how many cores there are -- ns
[ns,ms]=size(nms);

% Fill the string matrix indicating whether core already marked for excluding
strlive='-n';
strdie='-m';

% Initialize string matrix telling whether core to be kept or excluded
FV=repmat(strlive,ns,1);
% If cmask does not already exist, set to not mask any cores
if isempty(cmask) | ~(exist('cmask')==1);
   % no action needed;
else; % if cmask does exist, set the mask
   if any(~cmask); % if any cores previously marked for excluding
      nexclude = sum(~cmask);
      FV(~cmask,:)=repmat(strdie,nexclude,1);
   end;
end;



% Allow to build or revise the mask
dum=[];
for i=1:ms-4,
  dum=[dum,' '];
end
qt=['QUIT',dum];

k=0;
while k<ns+1,
  nmsn=[nms FV;qt,'  '];
  k=find(menudm2('Change core mask?',(cellstr(nmsn))'));
  if k<=ns,
     % Change the core mask for this core
     if cmask(k);
        cmask(k)=0;
        FV(k,:)=strdie;
     else;
        cmask(k)=1;
        FV(k,:)=strlive;
     end;
  end;
end;

% Make the mask logical
cmask=logical(cmask);
