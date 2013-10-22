function S=photo1
% photo1;
% photo1: surveys .jpg or .JPG files in current directory & allows interactive renaming
% Last Revised 1-13-01
%
% You have collected a tree-ring site and taken digital photos and want to QUICKLY process, and print out
% the photos.  The processing might include cropping, color adjustment, rotation, etc.
% You would like grayscales, several (say, 6) per page, labeled with the tree number for reference while 
% dating and interpreting the ring-width series.  You would also like the option of modifying the labels on the 
% "subimage" pages to include interesting information that might transpire (e.g., tacking on the dates).  You also want
% to save the full-size original .JPG images for possible later use individually (e.g., on a web page). Photo1,2,3,4
% are functions to simplify the task.  You could do these things with other software (e.g., Adobe Photodeluxe), but
% it would take much longer and your actions would not be easily repeatable.
%
%
%*** INPUT
%
% No arguments.
% User prompted to click on names of .jpg or .JPG files
%
%*** OUTPUT
%
% No arguments
% A file photolog.mat is created and stored in the current dir.  This file has a menu of filenames for later use.
% The original .JPG files are optionally deleted or renamed or left unchanged. Renaming allows to name photos with
%   sensible tree-name filenames.  For example, MHA01-1.JPG might be site MHA, tree 01, photo 1.
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE, except a subfunction
%
%*** TOOLBOXES NEEDED
% 
% Image processing (not in student group license)
% Last modified using version 6 of Matlab
%
%
%*** NOTES
%
% Usually Photo1, 2, 3 are run in sequence.  Photo4 and photo3 are run again if re-labeling of subimages desired.
%
%
%---Summary of steps for getting pages of multi-plotted figures
%
% Take photos of sample site with digigal camera
% Dump all the .JPG files from camera to empty directory
% Run photo1.m to rename or deleting .JPG files (also makes and saves photolog.mat)
% Run photo2.m to adjustphotos (e.g.,  , color change) --  optionally overwrites changed files
% Run photo2.m again to allow building of subimage figures -- also saves generation settings
% Run photo3.m to regenerate subimage page and print to file or printer
% Run photo4.m if desired to re-label the subimages (e.g., adding age info after trees dated)
% Run photo3.m again to regenerate page figure with revised labels and print figure
%
%
%--- Tips
% 
% Use the run of photo1 strictly to rename or delete files.
% Do all fine-tuning (rotation, color adjustment, cropping ...) in the first run of Photo2
% Make any pages of subimages in subsequent runs of photo2.m.  In these runs, the only thing you should
%    bother changing are the rgb2gray (if want grayscale subplots), and the scaling of the subplot images themselves
%    so that they fit well on the page and do not block the labels at the bottom of figure.  In these runs of
%    photo2,DO NOT OVERWRITE THE INDIVIDUAL .JPG FILES THAT HAD BEEN MODIFIED PREVIOUSLY.
% Use photo3.m to conveniently regenerate the subimage pages created interactively with Photo2.m. Photo3.m was 
%   written to avoid having to go back an duplictate all the interactive steps (photo selection, scaling, etc., 
%   used to originally make the subimage pages.
%
%--- Comment
%
% Some knowledge of graphics formats is necessary to make intelligent choices for image processing.  User should 
% read Matlab image-processing toolbox introductory material before trying these function.


% Get current directory
dir1 = cd;  % (1 x ?)s with directory
pf1 = [dir1 '\photolog.mat'];


%*** INITIAL FILE CHECK 
%
% Does photolog.mat exist in cd?  If yes, error. This should be initial look at photos.  
% If have a photolog.mat, will need to delete it.

close all;

%Check whether photolog.mat exists in cd
L = exist(pf1)==2;
if L; 
    error('photolog.mat already exists in current directory; cannot run photo1.m ');
else;  % photolog.mat not in cd
    [S1]=subfcn1(dir1,'JPG');
    [S2]=subfcn1(dir1,'jpg');
    S=[S1; S2];
end;
clear L kmen1 ;


% OK, if you are going on, that means you have run photo1.m before and stored S in photolog.mat, and you want to 
% accept the list of filenames and flags as store in S.


%***  BUILD MENU OF NAMES
nfile1 = size(S,1); % number of files
T=S;

% While loop over files
kwh1=1; 
while kwh1;
    T1=sort(T(1:nfile1));
    T=T1;
    T{nfile1+1}='Quit';
    for k=1:3;
        if ishandle(k);
            close(k);
        end;
    end;
        
    kmen1 = menu('Choose ',T);
    if kmen1==nfile1+1;
        kwh1=0;
    else;
        kwh1=1;
        f = T{kmen1};
        a=imread(f);
        aorig=a;
        figure(1);
        imshow(aorig);
        title([f ': original']);
        figure(2); 
        imshow(a);
        title([f ': current version']);
        
        % While loop to operate on files
        kwh2=1;
        while kwh2;
            kmen2=menu('Choose',...
                'Delete .jpg file from disk',...
                'Rename .jpg file',...
                'Return to file menu');
            switch kmen2;
            case 1; % delete file
                kquest = questdlg(['Are you sure you want to delete ' f ' from disk?']);
                if strcmp(kquest,'Yes');
                    
                    eval(['delete ' f ';']);
                    T{kmen1}=['Deleted: '  f];
                else;
                end;
           
                
            case 2; % rename file
                % Prompt for filename
                prompt={'Enter name:'};
                def={strtok(f,'.')};
                dlgTitle='New name for .jpg file (Do NOT include suffix)';
                lineNo=1;
                answer=inputdlg(prompt,dlgTitle,lineNo,def);
                if any(answer{1}=='.');
                    error('DUFUS, you entered a period in the name');
                end;
                
                fout=[upper(answer{1})  '.JPG'];
                
                
                % Write, checking that file with same name as writefile does not exist
                L=exist(fout,'file')==2;
                if L; % file of same name already exists
                    kquest=questdlg(['Sure you want to overwrite exising ' fout '?']);
                    if strcmp(kquest,'Yes');
                        eval(['!rename ' f  ' ' fout ';']);
                        T{kmen1}=[fout ' <-- ' f];
                    else; 
                    end;
                else; % file of same name does note yet exist
                    eval(['!rename ' f  ' ' fout ';']);
                    T{kmen1}=[fout ' <-- ' f];
                    
                    
                end;
               
            case 3; % Return fo file menu
                kwh2=0;
                
            end; % switch kmen2
            
            
            
        end; % while kwh2
    end;
end;





% SUBFUNCTIONS

function S=subfcn1(dir1,ext);
% Make cell vector of .jpg names
S=dirfls(dir1,ext,2); % String matrix of .jpg file names
S = cellstr(S); % cell vector of .jpg file names
n1 = size(S,1); % number of .jpg files in dirctory


function    [p1,scfact]=subfcn2(p1,scfact);
xlen = p1(3)*scfact;
ylen = p1(4)*scfact;
xadd = (xlen-p1(3))/2;
yadd = (ylen-p1(4))/2;xadd = (xlen-p1(3))/2;
yadd = (ylen-p1(4))/2;
p1(1)=p1(1)-xadd;
p1(3)=xlen;
p1(4)=ylen;








