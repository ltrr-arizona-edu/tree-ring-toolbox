function S=photo3
% photo3;
% photo3: automatic regeneration of page of subimages
% Last revised 1-17-01
%

%-- Load .mat file with needed setup info
close all;
[file1,path1]=uigetfile('page*.mat','infile of setup data as saved from photo2.m');
pf1=[path1 file1];
if exist(pf1,'file')~=2;
    error([pf1 ' does not exist ']);
end;
eval(['load ' pf1 ';']);

file1pref = strtok(file1,'.');


%--- Check contents


%--- Get arrangement for subimages


%--- Initialize subplot figure
% set position and paperposition
% turn off ticks
figure(1);

for n = 1:ncell;
    nm = phname{n};
    lab = phlabel{n};
    scfact=scalef(n);
    pcolor= phcolor{n};
    
    a = imread(nm);
    figure(2);
    imshow(a);
    figure(1);
    subplot(nrow,ncol,n);
    b=a;
    if ~isgray(b) & strcmp(pcolor,'Gray');
        b = rgb2gray(b);
    else;
    end;
    b=imresize(b,shrink,'bilinear');
    subimage(b);
    p1=get(gca,'Position');
    [p1,scfact]=subfcn2(p1,scfact);
    set(gca,'Position',p1);
    xlim=get(gca,'XLim');
    ylim=get(gca,'YLim');
    tt=text(xlim(1),ylim(2),lab,'FontSize',fslab,'VerticalAlignment','Top');
    
    
    set(gca,'XTick',[],'YTick',[]);
    
end;
set(gcf,'Position',pfig);
set(gcf,'PaperPosition',pappos);

%--  Print or export files

ptype={'jpeg','png','tiff','epsc2','meta','tiff','psc2'};
rlevel={'100','200','300','400','600','800','1200'}; % resolution in dpi
qual={'25','50','75','100'}; % quality level for jpeg

kmen1=menu('Choose output type',ptype);

str1 = ['-d' ptype{kmen1}] ;
str1b =' '; % tiff preview
str1a = ' ';
str1c=' '; % dpi
if kmen1==1; % jpeg
    kmen1a = menu('Choose JPEG quality (75 usual):',qual);
    str1a = [qual{kmen1a} ' '];
    
    
    
else;
end;
if kmen1==4; % epsc
    kquest=questdlg('Tiff preview?');
    switch kquest;
    case 'Yes';
        str1b='-tiff ';
    otherwise;
    end;
end;


str1 = [str1 str1a str1b];

% dpi
kmen2 = menu('Choose dpi (r level, 300 good',rlevel);
str2 = [' -r' rlevel{kmen2}]; 
    
    
% outfile
strdesc =[str1 str2];
if kmen1==6; % postscript, no file
    str3= ' ';
else;
    
    [file3,path3]=uiputfile(file1pref,[strdesc ' outfile']);
    pf3=[path3 file3];
    str3 = [' ' pf3];
end;


strall = ['print ' str1 str2 str3 ';'];

kq1=questdlg(['hgsave figure (will be as ' file1 ')?']);;
switch kq1;
case 'Yes';
    pf4=[path1 strtok(file1,'.')];
    eval(['hgsave ' pf4]);
otherwise;
end;

if kmen1~=7; % if not a psc2 
    kq2 = questdlg(['Print to file ' pf3 '?']);
    switch kq2;
    case 'Yes';
        eval(strall);
    otherwise;
    end;
else;
    kq2=questdlg(['Postscript print figure?']);
    switch kq2;
    case 'Yes';
        eval(strall);
    otherwise;
    end;
end;
    



function    [p1,scfact]=subfcn2(p1,scfact);
xlen = p1(3)*scfact;
ylen = p1(4)*scfact;
xadd = (xlen-p1(3))/2;
yadd = (ylen-p1(4))/2;xadd = (xlen-p1(3))/2;
yadd = (ylen-p1(4))/2;
p1(1)=p1(1)-xadd;
p1(3)=xlen;
p1(4)=ylen;


