function rwc2mat
% rwc2mat:  make an ascii time-series matrix from ring-width measurement structure
% rwc2mat;
% Last revise 10-22-01
%
%
%
%*** INPUT
%
% No args.  User prompted for files and settings
%
%*** OUTPUT 
%
% No args. User prompted for filename of output file.
%
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NOBE
%
%*** NOTES
%
% Function written to prepare V3-type input series for geosa1.m


%--- GET INPUT FILE

[file1,path1]=uigetfile('*.mat','Input .mat file (structure)');
pf1=[path1 file1];
eval(['s=load(''' pf1 ''');']);
if ~all([isfield(s,'vlist')  isfield(s,'XEL')  isfield(s,'XT')  isfield(s,'sinfo') isfield(s,'rwlset')]);
    error([pf1 ' does not contain XEL, XT , sinfo, rwlset and vlist']);
end
clear s;
% Load the ring-width data
eval(['load ' pf1 '  XT XEL sinfo rwlset vlist;']); 


%--- SELECT TYPE OF SERIES

kmen1=menu('Choose Series Type',...
    'Earlywood width',...
    'Latewood width',...
    'Total ring width');
if kmen1==1;
    xid = XEL.id;
    if isempty(xid{1});
        error(['No Earlywood or Latewood series in ' pf1]);
    else;
        jcol = 2;
        vbl='XEL';
    end;
elseif kmen1==2;
    xid = XEL.id;
    if isempty(xid{1});
        error(['No Earlywood or Latewood series in ' pf1]);
    else;
        jcol = 3;
        vbl='XEL';
    end;
else;
    xid = XT.id;
    if isempty(xid{1});
        error(['No Total Ring Width  series in ' pf1]);
    else;
        jcol = 2;
        vbl='XT';
    end;
end;




%--- BUILD MENU OF SERIES

eval(['x = ' vbl '.data']);
nser=length(xid);
h=repmat('-N',nser,1);
L=zeros(nser,1);
H=cellstr([char(xid) h]);
H{nser+1}='Finished';
kwh1=1;
while kwh1;
    kmen2=menu('Choose series',H);
    if kmen2<(nser+1);
        
        if strcmp(h(kmen2,:),'-N') & L(kmen2)==0;
            L(kmen2)=1;
            h(kmen2,:)='-Y';
        elseif strcmp(h(kmen2,:),'-Y') & L(kmen2)==1;
            L(kmen2)=0;
            h(kmen2,:)='-N';
        else;
        end;
        H=cellstr([char(xid) h]);
        H{nser+1}='Finished';
    else;
        kwh1=0;
    end;
    
end;


ipick = find(L);
npick = length(ipick);
if npick==0;
    error('None picked');
end;


% Find range of time series matrix
T=repmat(NaN,npick,2);
for n = 1:npick;
    ithis = ipick(n);
    d=x{ithis};
    yr = d(:,1);
    d = d(:,jcol);
    [d,yr]=subfcn01(d,yr);
    T(n,:)=[min(yr) max(yr)];
    plot(yr,d);
end;
yrgo = min(T(:,1));
yrsp = max(T(:,2));



%-- Allocate
yrX = (yrgo:yrsp)';
nsize = length(yrX);
X = repmat(NaN,nsize,npick);


%--- Fill Matrix;
for n = 1:npick;
    ithis = ipick(n);
    d=x{ithis};
    yr = d(:,1);
    d = d(:,jcol);
    [d,yr]=subfcn01(d,yr);
    i1 = min(yr) - yrgo +1;
    i2 = max(yr) - yrgo +1;
    X(i1:i2,n)=d;
end;



%---- SUBFUNCTION 01    
function [x,yrx]=subfcn01(x,yrx);
% Strip trailing and leading NaN
x=trailnan(x);
mx=length(x);
yrx=yrx(1:mx);
x=flipud(x);
yrx=flipud(yrx);
x=trailnan(x);
mx=length(x);
yrx=yrx(1:mx);
x=flipud(x);
yrx=flipud(yrx);
