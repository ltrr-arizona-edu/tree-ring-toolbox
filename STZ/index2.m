function index2
%
% DESCRIPTION : index2
% Converts tree indices in site chronologies, either standard or 
% residual
%
% INPUTS  :  NONE
% OUTPUTS :  NONE
%_________________________________________________

% Load the necessary .mat file(s)
% Get the .mat filename interactively
flmat=uigetfile('*.mat','.MAT filename ?');
eval(['load ',flmat]);

% Check if the mat file exists in the workspace
if ~(exist('IT') & exist('ET')),
  udisp('Wrong file. Please enter correct filename');
  flmat=uigetfile('*.mat','.MAT filename ?');
  eval(['load ',flmat]);
end

% Save the original variables in a temporary file
save tempo.mat;
