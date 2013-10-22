function coremask
%
% DESCRIPTION : coremask
% This function builds the mask for user to specify whether or 
% not to include core indices in the chronology. User will base 
% decision on statistics and other information gained so far. Only
% the unmasked series will be used in subsequent analyses.
%
% INPUTS  :  NONE
%
% OUTPUTS :  NONE
%_________________________________________________________________



% Prompt for the .mat file name
% Load the necessary .mat file(s)
flmat=uigetfile('*.mat','.MAT filename ?');
eval(['load ',flmat]);

% Check if nms matrix exists
if ~exist('nms'),
  jh=jdisp('Core ID matrix does not exist! Wrong filename.');
  pause;
  close(jh);
  return;
end

% Save the original variables in a temporary file
save tempo.mat;

% Initialize the mask identifier
[ns,ms]=size(nms);
fv=[];
for i=1:ns,
  fv=[fv;'-n'];
end

if exist('cmask'),
  mu1=usinp('core mask (cmask) already exists! Overwrite ?');
  if ~mu1,
    return;
  end
end

% Build the mask
dum=[];
for i=1:ms-4,
  dum=[dum,' '];
end
qt=['QUIT',dum];
k=0;
cmask=ones(ns,1);
while k<ns+1,
  nmsn=[nms,fv;qt,'  '];
  k=svmenu('Mask Core ID?',nmsn);
  if k<=ns,
    cmask(k)=0;
    fv(k,:)='-m';
  end
end

% Display the mask vector
disp(cmask);

% Prompt for the file name to store cmask
nsv=shlmnu('Select','Add to file?','Store in new file?','QUIT');
if nsv==1,
  eval(['save ','temp.mat',' cmask']);
  clear all;
  load tempo;
  load temp;
  eval(['save ',flmat]);
elseif nsv==2,
  ofmat=uiputfile('*.mat','.MAT filename ?');
  eval(['save ',ofmat,' cmask']);
end


% End of file
