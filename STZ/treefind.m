function [treenms,itree]=treefind(corenms);
% treefind:   Extract tree names from character matrix of core names
% [treenms,itree]=treefind(corenms);
% Last revised 2-22-02
%
% Extract tree names from character matrix of core names. Called by lfcheck.m to help in grouping matrix
% of core indices into tree indices.
%
%*** INPUT
%
% corenms (? x ?)s   core names, in specified convention of naming (see notes)
%
%
%*** OUTPUT
%
% treenms{}s     coll-cell with sting names of trees
% itree (? x ?)i  index to rows of corenms indicating which tree (in treenms) the core belongs to
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED  -- NONE
%
%*** NOTES
%
% Assumed format of core names in corenms:
%   First 3 chars (disregarding leading blanks) are the site code, either letters or numbers
%   Next 1-3 chars are tree number; must be numeric, in range 1-999
%   Next char is the core id; mist be a letter, and is required
%   Next (optional) char is a numeric within-core sequence number, numbered from inside out


C=corenms;
[mC,nC]=size(C);
if ~ischar(C);
    error('C must be character matrix');
end;
T=repmat(blanks(3),mC,1); % to store tree  number


% Conver to upper case
C=upper(C);

%--- LEFT JUSTIFY char matrix
C=strjust(C,'left');

%-- Strip off site codes
S=C(:,1:3); %  store site codes
C1=C;
C(:,1:3)=[];

%-- Form 3-char all-number tree-number string matrux
for n=1:mC; % loop over core names
    c=C(n,:);
    i1 = find(isletter(c));
    if isempty(i1);
        error(['No core letter for core name ' corenms(n,:)]);
    end;
    if (i1==1);
        error(['Fourth character of core name ' corenms(n,:) ' not numeric']);
    end;
    c1 = c(1:(i1-1));  % pull numeric part
    x1=str2num(c1); % tree number as numeric
    if x1<=0;
        error(['Tree number of core  ' corenms(n,:) ' not greater than 0']);
    end;
    if x1<10;
        c1=['00' num2str(x1)];
    elseif x1<100;
        c1=['0' num2str(x1)];
    elseif x1<1000;
        c1=num2str(x1);
    else;
        error(['Tree number of core  ' corenms(n,:) ' greater than 999']);
    end;
    if length(c1) ~=3;
        error(['Length of ' c1 ' not 3']);
    end;
    T(n,:)=c1;
end;

% Trim leading digit off tree numbers if all less than 100
if  all(T(:,1)=='0');
    T(:,1)=[];
end;


% Merge site code and tree number
T=[S T];


% Assign cores to trees

kwh1=1;
j0=(1:mC)';
ntree=0;
igo=1;
ithis=igo;
J=(2:mC)';
itree=repmat(NaN,mC,1);
while kwh1;
    ntree=ntree+1; % current # trees identified
    w=T(ithis,:);
    W=cellstr(T(J,:));
    I=strmatch(w,W,'exact');
    itree(ithis)=igo;
    treenms{igo}=w;
    if ~isempty(I);
        itree(J(I))=igo;
        J(I)=[];
    else;
    end;
    if ~isempty(J);
        if length(J)==1;
            itree(J)=igo+1;
            treenms{igo+1}=T(J,:);
            kwh1=0;
        else;
            ithis=J(1);
            J(1)=[];
            igo=igo+1;
        end;
    else;
        kwh1=0;
    end;
end;

treenms=treenms'; % convert to col-cell