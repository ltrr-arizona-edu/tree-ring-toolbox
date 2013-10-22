function S=jpgsize

%*** NOTES
%
% cd to directory with photo files before running
%
% 

kmen1=menu('Choose Target Width (pixels)','600','450');
if kmen1==1;
    wspec=600;
elseif kmen1==2;
    wspec=450;
end;



% Get current directory
dir1 = cd;  % (1 x ?)s with directory
pf1 = [dir1 '\photolog.mat'];


%*** INITIAL FILE CHECK 
%
% Does photolog.mat exist in cd?  If not, create photolog.mat and store a cell vector of filenames S in it; then quit.
% If yes, does photolog.mat contain S?  If no, error message.  If yes, prompt for whether you want to use the old cell vector of 
% names or build and store a new one.  The only way you continue in this function is if you accept the existing S.

close all;

nwind=4;
nsub=6;
jsub=0;

%Check whether photolog.mat exists in cd
L = exist(pf1)==2;
if L; 
    load photolog;
    if exist('S','var')==1; 
        kmen1=menu('Choose',...
            'Use the existing cell vector of filenames S',...
            'Re-make the cell vector of filenames S');
        switch kmen1;
        case 1; % use existing cell vector or names
        case 2; % Make new cell vector S, store it, and exit program
            [S]=subfcn1(dir1,'JPG');
            save photolog S;
            %return;
        end; % switch kmen1
    else; % 
        disp([pf1 ' exists in cd, but does not contain S']);
        error('See messsage above');
    end;
else;  % photolog.mat not in cd
    [S]=subfcn1(dir1,'JPG');
    save photolog S;
    %return;
end;
clear L kmen1 ;



% OK, if you are going on, that means you have run photo1.m before and stored S in photolog.mat, and you want to 
% accept the list of filenames and flags as store in S.


%***  BUILD MENU OF NAMES
nfile1 = size(S,1); % number of files
T1=repmat(NaN,nfile1,2); % to hold original size
T2=repmat(NaN,nfile1,2); % to hold Web size

for n = 1:nfile1; % loop over files
    s = S{n};
    A=imread(s);
    s3  = size(A); % 1 x 3, with width as dimension 2
    s3(3)=[];
    T1(n,:)=s3;
    rat1 = wspec/s3(2); % ratio of desired width to original width
    T2(n,:)=s3*rat1;
end;
T2=round(T2);
T3 = [char(S) repmat(blanks(2),nfile1,1) (num2str((T1))) repmat(blanks(2),nfile1,1)  (num2str((T2)))];

uiwait(msgbox('Copy, paste, print, etc. as needed form the command window','Message','modal'));

clc;
disp(T3);





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








