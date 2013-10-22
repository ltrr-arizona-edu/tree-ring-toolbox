function pre=dirfls4(path1,suff)
% dirfls4: list files with specified suffix(es) are in a directory
% pre=dirfls4(path1,suff);
% Last revised 10-10-01
%
% List files with specified suffix(es) are in a directory.  Example of use: find all .eww and lww files
% in c:\work5\.
% 
%
%*** INPUT
%
% path1 (1 x ?)s  directory <'c:\work5\'>
% suff {} cell array of suffixes <{'eww','lww'}
%
%
%*** OUTPUT 
%
% pre. structure with fields x1, x2, ... 
%   Each field has a cell array of files with the suffix in suff{1}, suff{2}, etc,
%   or is returned empty (see notes)
%   
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLE -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% pre.x1{} on output holds the prefixes of all files with suffix suff{1}.  pre.x2 ... with suffix suff{2}, etc.
% pre.x1 is returned empty if no files with suffix suff{1}, etc.



nsuff=size(suff,2); % number of suffixes
pre=cell(nsuff,1);

D=dir(path1);  % puts file info into structure D
nfls1 = size(D,1);  % number of files, including directories, in current dir
C=struct2cell(D); % convert structure to cell. C(1,:) holds file names
c = C(1,:); % file names, in cell array


for j = 1:nsuff; % loop over suffixes
    clear prefx;
    
    suffx=suff{j};
    lensuff=length(suffx);
    ncount=0; 
    for n = 1:nfls1; % Loop over filenames
        d1 = c{n}; %  a single filename, as char -- could also be a directory
        d1 = strtok(d1); % truncate trailing blanks
        i1 = findstr(d1,lower(suffx));
        i2 = findstr(d1,upper(suffx));
        if ~(isempty(i1) & isempty(i2));
            if isempty(i1);
                i1=i2;
            else
                i1=i1;
            end
            i1=i1(length(i1)); % to make sure pointing to start of last occurrence of suff
            % Make sure suffix contains exactly the right number of characters
            if length(d1)==i1+lensuff-1;
                ncount=ncount+1;
                d1 = strtok(d1,'.'); % drop off the suffix
                prefx{ncount}=d1;
                
            else
            end
        end
        
    end % for over filenames
    if exist('prefx')==1;
        eval(['pre.x' int2str(j) '=prefx;']);
    else;
        eval(['pre.x' int2str(j) '=[];']);
    end;
end; % for over suffixes

if isempty(pre);
    error(['No files found with any of suffixes : ' suff]);
end


