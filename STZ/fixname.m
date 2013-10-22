function fixname
% fixname: fix unconventional core names in a .mat storage matrix
% fixname;
% Last revised 12-18-00
%
% Fixes invalid core ids (names) so that crvfit1a, crvfit1b, etc. work properly. Needed because of irregular
% naming conventions.  For example, some latewood cores have ids ending in XL (e.g., JAC05AXL.RW).  
%
%*** INPUT
%
% No args.
% User prompted to click on .mat storage file holding ringwidths in column-vector format in X and
%    core names in nms
%
%*** OUTPUT
%
% No output args. 
% User prompted whether to store revised names in the input .mat file, overwriting the 
%   old nms.  Option to store in new .mat file.
%
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% A valid names should the following form:
%   3 leading letters (site code)
%   A number (tree number)
%   A Letter (core designator)
%   <a number> portion number within core (numbered from inside out)

% GET THE FILE WITH RINGWIDTHS
[file1,path1]=uigetfile('*.mat','Infile with ringwidths and nms');
pf1=[path1 file1];
eval(['load ' pf1 ' nms;']);

nmsold=nms; % store original names
[m1,n1]=size(nms);
a = nms(:,1:3); % first 3 chars
L=isletter(a);
if ~all(all(L'));
    disp(nms);
    error('Look, first 3 chars not all letters for all series');
end;

b=nms(:,4:n1); % remainder of id, beginning with tree number

% Loop over series
for n = 1:m1;
    d = b(n,:);
    L = isletter(d);
    i=find(L);
    if min(i)<2;
        disp(d);
        error('Look, no tree number');
    else;
        dnum = d(1:(min(i)-1));
        j = str2num(dnum);
        if isempty(j);
            disp(d);
            error('Look, tree number is not all number');
        end;
    end;
    
        
        nletter= sum(L);
    if nletter==0;
        disp(d);
        error('Look, this part of name has no letter following tree number');
    elseif nletter>1;
        if ~all(diff(i)==1);
            disp(d);
            error('Look, letters after first 3  chars not contiguous');
        end;
        if nletter==3;
            e = upper(d(i(2):i(3)));
            if strcmp(e,'XL') | strcmp(e,'XE');
                e = blanks(2);
                d(i(2):i(3)) =e;
                b(n,:)=d;
            else;
                disp(d);
                error('Look, name has 3 letters, but end is not XL or XE');
            end;
        else;
            disp(d);
            error('Look, name has 2 letters here');
        end;
    else;
        % one trailing letter, OK
    end;
    
end;
% 
% 
nms=[a b];

if strcmp(nms,nmsold);
    disp('No change to names needed');
else;
    
    B=[nmsold  repmat(blanks(3),m1,1)  nms];
    disp('Old and revised names');
    
    % 
    disp(B);
    
    eval(['save ' pf1 ' nms -append;']);
end;

