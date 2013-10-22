function corev

% Intialize the ui
uicorev; % calling with no args initializes

clear all

% Prompt for the .mat file name

% Get the .mat filename interactively; load file
flmat=uigetfile('*.mat','.MAT filename ?');
flold=flmat; % Will need this file name later
eval(['load ',flmat]);


if ~exist('WE') | ~exist('WI') | ~exist('rbtE') | ~exist('rbtI'),
	error('WE, WI, rbtE, or rbtI not in .mat file')
end


global WE WI rbtE rbtI


epsv
