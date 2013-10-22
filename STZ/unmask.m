function unmask
% unmask: remove masking or blocking of core series
% unmask;
% Last revised 4/4/01
%
% You had masked out cores from the chronology using crvfit.  Or you had blocked out
% segments of some cores.  Now you want to reverse the masking or blocking.  unmask.m allows you 
% remove the masks or blocks globally (all series), or individually.  You MUST re-run the 
% sequence corei, treei, etc., after running unmask.
%
%*** INPUT
%
% No args.
% User prompted for name of input .mat chronolology storage file
%
%*** OUTPUT
%
% No args
% User prompted for name of output .mat chronology storage file.  Default is
% to replace the input file
%
%*** REFERENCE -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% Recall that:
%   cmask is a cv with 1  or zero for each series.  0 means core masked
%   S(i,13) is 0 or 1 depending on whether blocks imbedded in series i
%   cmaskwhy{i,3}
%   blockwhy{i x 3} cell blocking information, series i
%       {i,1} 2-col matrix with first, last year of block
%       {i,2} char vector letter code of reason for each block, in order as in {i,1}
%           H growth anomalously high
%           L growth anomalously low
%           C correlation with other series weak (e.g., by COFECHA)
%           
%


%--- GET INPUT DATA

[file1,path1]=uigetfile('*.mat','Input .mat chronology storage file');
pf1=[path1 file1];

% Load the .mat storage file
eval(['load ' pf1]);
if ~exist('nms','var')==1 | ~exist('S','var')==1 | ~exist('yrs','var')==1;
    error(['nms, S, yrs  not in ' pf1]); 
end;
ncores = size(nms,1);

%--- UNMASKING

kmen1=menu('Choose one',...
    'Remove masks individually',...
    'Remove masks globally',...
    'No unmasking, move on to unblocking',...
    'Abort and bail out');

if kmen1==1 | kmen1==2; % unmasking
    if ~exist('cmask','var')==1;
        error(['cmask not in ' pf1]); 
    end;
    Lmasked=cmask==0; % logical index to currently masked series
    sum1=sum(Lmasked);
    if sum1==0;
        error('No series currently masked -- nothing to unmask');
    end;
end;

vnms = nms;
if kmen1==1; % individual unmasking
    
    d = repmat(blanks(2),ncores,1);
    d(Lmasked,1)='-';
    dnew='  ';
    d(Lmasked,2)='X';
    vnms=[vnms d];
    cname=cellstr(vnms);
    cname{ncores+1}='No more unmasking';
    
    kwh1=1;
    while kwh1==1;
        kmen2=menu('Choose',cname);
        if kmen2==(ncores+1);
            kwh1=0;
        else;
            if cmask(kmen2)==1; % no action, not masked
            else;
                cmask(kmen2)=1;
                cname{kmen2}=[nms(kmen2,:) dnew];
            end;
        end;
    end;
elseif kmen1==2; % global unmasking
    cmask=ones(ncores,1);
elseif kmen1==3; % move on to unblocking
else; % abort
    return;
end;


pf2def = pf1;
kmen3=menu('Choose',...
    ['Save revised mask/block information in input file ' upper(pf2def)],...
    'Save revised mask/block and all other data in a new .mat file',...
    'Abort');
if kmen3==1;
    eval(['save ' pf2def ' cmask -append;']);
elseif kmen3==2;
    pf2=uiputfile('*.mat','Save revised masking/blocking information here');
    eval(['save ' pf2  ' cmask;']);
else;
end;



disp('Starting to unblock -- a future implementation');



    
                
        
    
    
    





    
    