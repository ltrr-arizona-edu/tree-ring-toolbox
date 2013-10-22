function clisplc2
% clisplc2: following clisplc1.m, compare the augmented with the old series; optional update
% clisplc1;
% Last revised 10-6-99
%
%*** INPUT
%
% You are prompted to click on file names:
%   input old series
%   input updated series
%   output updated series
%
%*** OUTPUT
%
% If you approve, the updated series is accepted, and is optionally saved with a new 
% filename, and optionally the old input series and the recent-data file are
% deleted.  Expected naming sequence, assuming root name of tucs:
%
%  tucs.mat --- the input old series
%  tucs99.mat -- the recent data
%  tucsnew.mat -- the spliced data (from running clisplc1.m
%
% Typically, if all looks OK, tucs98.mat would be deleted, tucsnew.mat would be 
% saved as tucs.mat, and then tucsnew.mat would be deleted.
%
%*** REFERENCES --- NONE
%*** UW FUNCTIONS CALLED -- none
%*** TOOLBOXES NEEDED -- none
%
%*** NOTES


%--- Click on the old input file
[file1,path1]=uigetfile('*.mat','Old (before updating) .mat file');
pf1=[path1 file1];

%--- Load old data
eval(['load ' pf1 ' Z;']);
Zold=Z;

%--- Load tentative updated data
fname = strtok(file1,'.');
[file2,path2]=uigetfile([path1 fname 'new.mat'],'Infile with updated data');
pf2=[path2 file2];
eval(['load ' pf2 ' Z;']);


%--- Compare each month's data
for i =1:12;
   hp1=plot(Zold(:,1),Zold(:,i+1));
   hold on;
   hp2=plot(Zold(:,1),Zold(:,i+1),'*');
   hp3=plot(Z(:,1),Z(:,i+1));
   hp4=plot(Z(:,1),Z(:,i+1),'o');
   set(hp1,'LineWidth',3,'Color',[.5 .5 .5]);
   set(hp2,'Color',[.5 .5 .5]);
   set(hp3,'Color',[1 0 0]);
   set(hp4,'Color',[1 0 0]);
   hold off;
   str1=[fname ': Month ' int2str(i)];
   title(str1);
   grid;
   zoom xon;
   pause;
   
end;

%--- Optionally overwrite old file
kq1=questdlg('Accept revised data?');
if strcmp(kq1,'Yes')
   eval(['save ' pf1 ' Z;']);
   eval(['delete ' pf2 ';']);
   f98=[fname '99.mat'];
   pf3=[path2 f98];
   eval(['delete ' pf3]);
end;




   





%   