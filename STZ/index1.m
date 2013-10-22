function index1
%
% DESCRIPTION : index1
% Converts core indices into tree indices
%
% INPUTS  :  NONE
% OUTPUTS :  NONE
%_________________________________________________



% Load the necessary .mat file(s)
% Get the .mat filename interactively
flmat=uigetfile('*.mat','.MAT filename ?');
eval(['load ',flmat]);

% Check if the mat file exists in the workspace
if ~(exist('EX') & exist('IX')),
  udisp('Wrong file. Please enter correct filename');
  flmat=uigetfile('*.mat','.MAT filename ?');
  eval(['load ',flmat]);
end

% Save the original variables in a temporary file
save tempo.mat;

% Call treeno for tree names 
[Inms,tnms,n]=treenum(nms,cmask);

% Call meantree to compute tree indices
[IT,ITn,Tnms,ITyrs] = meantree(IX,nms,yrs,cmask);
[ET,ETn,Tnms,ETyrs] = meantree(EX,nms,yrs,cmask);

% Save the vectors in a .mat file
nsv=shlmnu('Save variables?','Add to file?','Store in new file?','QUIT');
if nsv==1,
  eval(['save ','temp.mat',' IT ET ITn ETn ITyrs ETyrs Tnms']);
  clear all;
  load tempo;
  load temp;
  eval(['save ',flmat]);
elseif nsv==2,
  ofmat=uiputfile('*.mat','.MAT filename ?');
  eval(['save ',ofmat,' IT ET ITn ETn ITyrs  ETyrs Tnms']);
end


% End of file
