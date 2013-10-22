function treei
%
% DESCRIPTION : treei
% Converts core indices into tree indices
%
% INPUTS  :  NONE
% OUTPUTS :  NONE

%*** UW FUNCTIONS CALLED
%
% coremask
% treenum
% meantree
% stats2


clear all

% Prompt for the .mat file name

% Get the .mat filename interactively; load file
[flmat,path1]=uigetfile('*.mat','Input original .MAT filename ?');
pf1=[path1 flmat];
flold=flmat; % Will need this file name later
eval(['load ',pf1]);
clear flmat;

% Check if the mat file exists in the workspace
if ~(exist('EX') | ~exist('IX')),
  error('This .mat file does not have EX and IX')
end


% Build a mask to omitt some cores if desired -- 0 means omitt
if ~(exist('cmask')==1), 
   cmask=[]; 
   cmask=logical(cmask);
else;
   cmask=logical(cmask); 
end
cmask = coremask(nms,cmask)



% Call treenum for tree names 
[Inms,tnms,n]=treenum(nms,cmask);

% Call meantree to compute tree indices
[IT,ITn,Tnms,ITyrs] = meantree(IX,nms,yrs,cmask);
[ET,ETn,Tnms,ETyrs] = meantree(EX,nms,yrs,cmask);


% Statistics on tree indices
[rbtI, DI, WI]=stats2(IT,Tnms,ITyrs);
[rbtE, DE, WE]=stats2(ET,Tnms,ETyrs);


% Save the vectors in a .mat file
newvars1=' cmask IT ET ITn ETn ITyrs ETyrs Tnms '; % new variables
newvars2=' rbtI rbtE DI DE WI WE '; % more new variables
newvars=[newvars1 newvars2];
nsv=menu('Save variables?','Add to Original .mat file?',...
   'Store new variables in new .mat file?','QUIT');
if nsv==1,
  eval(['save ' pf1 newvars ' -append']);
elseif nsv==2,
  [ofmat,path2]=uiputfile('*.mat','.MAT file to store new vars in: ');
  pf2=[path2 ofmat];
  eval(['save ' pf2 newvars]);
end

% End of file
