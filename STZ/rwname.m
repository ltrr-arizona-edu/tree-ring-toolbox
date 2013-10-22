function rwname
% rwname:  re-name .rw,.eww, or .lww files
% rwname;
% Last revised 4-4-01
%
% You have, say, latewood files with irregular names, such as dcy9axe.eww.  You prefer the
% more consistent dcy09a.eww, which has
%   3-letter site code
%   tree number, preferably with "0" before single digits to facilitate file listing
%   core letter
%   optional core-segment number:  e.g., dcy09a1.eww
%
%*** INPUT
%
% No args
% User assembles files of desired suffix (e.g., *.eww) in an empty directory.  rwname.m
% operates on these files & writes re-named files to another directory
%
%*** OUTPUT 
%
% The files with revised names
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES



% Prompt for file type
kmen1=menu('Choose',...
    '.rw files',...
    '.eww files',...
    '.lww files');
switch kmen1;
case 1;
    suffx='rw';
case 2;
    suffx='eww';
case 3;
    suffx='lww';
otherwise;
end;


% Prompt for output directory
[file1,path1] = uigetfile('*.*','Click on any file to get target directory for output (renamed) files');
% path1 is output dir


% PROMPT FOR OPTIONAL CHANGE IN 3-LETTER CODES
kquest1=questdlg('Change any 3-letter codes');
switch(kquest1);
case 'Yes';
    kcode=1;
otherwise;
    kcode=0;
end;

if kcode==1;
    prompt={'Enter the original code:','Enter the revised code:'};
    def={'XXX','xxx'};
    dlgTitle='Option code change (e.g. JMD to jmr)';
    lineNo=1;
    answer=inputdlg(prompt,dlgTitle,lineNo,def);
    codeold = answer{1};
    codenew = lower(answer{2});
end;



%--- SURVEY DIRECTORY AND STORE RAW FILENAMES IN A STRUCTURE

eval(['p = dir(''*.' suffx ''');']);
% p.names now has the file names

nfiles=size(p,1);


%-- Rename

for n = 1: nfiles;
    d=p(n).name;
    d1=strtok(d,'.'); % file prefix
    
    % First 3 chars should be letters
    L1 = isletter(d1);
    if not(all(L1(1:3)) & not(L1(4)));
        error([d1 ' does not have 3 letters follwed by number']);
    else;
        scode=lower(d1(1:3));  % lower case of site code
        
        if kcode==1;
            if strcmp(scode,codeold);
                scode=codenew;
            end;
        end;
        
        d1(1:3)=[]; % strip off 3 letter code
        L1=isletter(d1);
        diff1 = diff(L1);
        i1 = find(diff1==1)+1; % position of next letter ( core letter);
        n1 = i1-1; % number of characters in tree number
        dtree = d1(1:n1);
        if n1 ==1;
            stree = ['0' (dtree)];
        else;
            stree = (dtree);
        end;
        d1(1:n1)=[] ; % strip off tree number
        if isempty(d1);
            slett = 'a';
        else;
            slett = lower(d1(1));
            d1(1)=[]; % strip off core letter;
            L1 = isletter(d1);
            if all(L1); % no within core seg number
                seg='No';
            else;
                nzero = sum(~isletter(d1)); % number of numbers after core letter
                if nzero==1 & ~isletter(d1(length(d1)));
                    ssegment = (d1(length(d1)));
                    seg='Yes';
                    
                else;
                    error(['undecipherable form ' d]);
                end;
            end;
            
            if strcmp(seg,'Yes');
                fn = [scode stree slett ssegment '.' suffx];
            else;
                fn = [scode stree slett '.' suffx];
            end;
            
            pf1 = [path1 fn];
           
        end;
        
        disp([d ' to ' pf1]);
        
      
    end;
    copyfile(d,pf1);
end;
    
    
    
    
   
    
    