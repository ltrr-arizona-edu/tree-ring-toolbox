function [z,id] = get1ts

% USAGE : [z,id] = get1ts
% This function gets one time series out of a .mat file
%
% INPUTS	: NONE
%
% OUTPUTS	: z (2 x n) - desired time series
%		  id (String) - Core ID
%
%_______________________________________________________________



% Prompt for the .mat file name

% Get the .mat filename interactively
flmat=uigetfile('*.mat','.MAT filename ?');

% Check if the mat file exists in the workspace
if ~(exist('X') & exist('IX')),
  eval(['load ',flmat]);
end

[ns,dum1]=size(nms);
nsmt=[];
for i=1:ns,
 nsmt=str2mat(nsmt,[int2str(i),'-']);
end
nsmt(1,:)=[];

nsms=[nsmt nms];
% Prompt for the core ID sequence
scid=svmenu('Core ID # ?',nsms);

% Specify the desired time series matrix
sx=shlmnu('Pick Series','X','IX','EX','Other');
yrv=(yrs(scid,1):yrs(scid,2))';
if sx==1,
  xv=X(yrs(scid,3):yrs(scid,3)+yrs(scid,2)-yrs(scid,1));
elseif sx==2,
  xv=IX(yrs(scid,3):yrs(scid,3)+yrs(scid,2)-yrs(scid,1));
elseif sx==3,
  xv=EX(yrs(scid,3):yrs(scid,3)+yrs(scid,2)-yrs(scid,1)); 
else
  jh=jdisp('Please return to command window for input');
  pause(1);
  close(jh);
  snm=input('Please enter the matrix name = ','s');
  if exist(snm),
    eval(['xv=',snm,'(yrs(scid,3):yrs(scid,3)+yrs(scid,2)-yrs(scid,1));']);
  else
    jh=jdisp('Specified matrix does not exist. Please rerun the program');
    pause(2);
    close(jh);
    return;
  end
end
z=[yrv,xv];
id=nms(scid,:);

 % Plot the rw data with scid core ID
 hf0=figure('Color','w','Units','normal');
 plot(yrv,xv,'b');title(['RingWidth ',nms(scid,:)]);
 xlabel('Year');ylabel('RW, (0.01 mm)');


% End of file
