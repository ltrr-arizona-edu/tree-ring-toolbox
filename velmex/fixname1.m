function nms=fixname1(nms)
% fixname1: fix unconventional core names in a .mat storage matrix
% fixname1;
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



nmsold=nms; % store original names
[m1,n1]=size(nms);
P=repmat(NaN,m1,1); % to store numeric tree numbers
B=cell(m1,1); % to store part of names after tree number
a = nms(:,1:3); % first 3 chars
% L=isletter(a);
% if ~all(all(L'));
%     disp(nms);
%     error('Look, first 3 chars not all letters for all series');
% end;

b=nms(:,4:n1); % remainder of id, beginning with tree number

%--- Check for old core-numbering convention that used 1,2,3, ... instead of a, b, c, ... for core labels
L=isletter(b);
if ~any(any(L));
    u = {b};
    u=deblank(u);
    u1=char(u);
    u2=u1(:,end);
    L1=u2=='1' | u2=='2' | u2=='3' | u2=='4' | u2=='5';
    if all(L1);
        koldconv=1;
        h = {'a','b','c','d','e'};
        for n =1:m1;
            bthis=b(n,:);
            i1=find(~isspace(bthis));
            j = str2num(bthis(max(i1)));
            bthis(max(i1))=h{j};
            b(n,:)=bthis;
        end;
    else;
        koldconv=0;
    end;
end;


% Loop over series
for n = 1:m1;
    d = b(n,:);
    L = isletter(d);
    i=find(L);
    if isempty(i);
        disp(nms(n,:));
        error('Look, no letter following tree number');
    end;
    if min(i)<2;
        disp(nms(n,:));
        error('Look, char 4 is not a number');
    else;
        f=upper(d(min(i):end)); % remainder of name after tree number
        dnum = d(1:(min(i)-1)); % tree number
        j = str2num(dnum);
        P(n)=j;
        
        L=isletter(f);
        i=find(L);
        nletter= sum(L);
        if nletter==0;
            disp(nms(n,:));
            error('Look,  no letter following tree number');
        elseif nletter>1;
            if ~all(diff(i)==1);
                disp(nms(n,:));
                error('Look, letters after tree number not contiguous');
            end;
            if nletter==3;
                e = upper(f(i(2):i(3)));
                if strcmp(e,'XL') | strcmp(e,'XE');
                    f(i(2):i(3)) =[];
                else;
                    disp(nms(n,:));
                    error('Look, name has 3 letters, but end is not XL or XE');
                end;
            else;
                disp(nms(n,:));
                error('Look, number of letters after tree number either 2 or greater than 3');
            end;
        else;
            % one trailing letter, OK
        end;
    end;
    B{n}=f;
end;
B=char(B);


% Adjust strings of tree numbers to uniform number of chars
Pstr=cell(m1,1);
nummax = max(P);
for n=1:m1;
    p=P(n);
    if nummax<10;
        Pstr{n}=num2str(p);
    elseif nummax<100;
        if p<10;
            Pstr{n}=['0' num2str(p)];
        else;
            Pstr{n}=num2str(p);
        end;
    else;
        if p<10;
            Pstr{n}=['00' num2str(p)];
        elseif p<100;
            Pstr{n}=['0' num2str(p)];
        else;
            Pstr{n}=num2str(p);
        end;
    end;
end;
Pstr=char(Pstr);


nms=upper([a Pstr B]);
