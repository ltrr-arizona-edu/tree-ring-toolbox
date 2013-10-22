function S=photo4
% photo4;
% photo4: revise labeling of subimages on figure produced by photo2.m
% Last revised 1-16-01
%
% photo2.m optionally built a figure of subimages, and interactively allowed labeling of 
% the images.  photo3.m allows user to modify the labels on the subimages.  photo3.m assumes the
% figure was save with hgsave from photo2.m.   photo3.m loads that file with hgload, finds the handles to the text 
% objects, and interactively allows those objects to be modified.
                        

%-- Load the saved figure file


close all;
[file1,path1]=uigetfile('*.fig','.fig infile of previously labeled subplots');
pf1=[path1 file1];
if exist(pf1,'file')~=2;
    error([pf1 ' does not exist ']);
end;
h=hgload(pf1);
pf2=strtok(pf1,'.');
eval(['load ' pf2 ' phlabel fslab ;']);



% Prompt whether want to overwrite font 
kfont = questdlg('Specify new uniform font size?');
switch kfont;
case 'Yes';
    prompt={'Enter font size:'};
    def={'10'};
    dlgTitle='Desired font size';
    lineNo=1;
    answer=inputdlg(prompt,dlgTitle,lineNo,def);
    font1= str2num(answer{1});
case 'No';
    font1=NaN;
case 'Cancel';
    font1=NaN;
otherwise;
end;

    

% Find the children of the figure
hc=get(gcf,'Children');
hc = findobj(hc,'Type','Text');
fontold = get(hc(1),'FontSize');
nfig = length(hc);

% Store the text strings in cell
a = cell(nfig,1);

for n = 1:nfig;
    str1 = get(hc(n),'String');
    a{n}=str1;
end;
T=a;

T{nfig+1}='Quit';

% While loop over subfigures


kwh1=1;
while kwh1;
    kmen1 = menu('Choose one',T);
    if kmen1==(nfig+1);
        kwh1=0;
    else; 
        kwh1=1;
        str1=T{kmen1};
        hthis = hc(kmen1);
        deftxt =str1;
        
        switch kfont;
        case 'Yes';
            prompt={'Enter revised text:','Enter the fontsize:'};
            def={deftxt,num2str(font1)};
            dlgTitle='Modify text';
            lineNo=1;
            answer=inputdlg(prompt,dlgTitle,lineNo,def);
            newtext = char(answer{1});
            newfont = str2num(answer{2});
            T{kmen1}= newtext;
        otherwise;
            prompt={'Enter revised text:'};
            def={deftxt};
            dlgTitle='Modify text';
            lineNo=1;
            answer=inputdlg(prompt,dlgTitle,lineNo,def);
            newtext = char(answer{1});
            newfont = fontold;
            T{kmen1}=newtext;
            
        end; % switch kfont
        
        
        % Modify
        set(hthis,'String',newtext,'FontSize',newfont);
        phlabel{nfig-kmen1+1}=newtext;
              
    end;
end;

% Prompt to overwrite .fig file or abort
kmen2=menu('Choose',...
    ['Overwrite ' pf1 ' with revised version of figure'],...
    'Abort -- no overwrite');
switch kmen2;
case 1;
    eval(['hgsave ' pf1]);
    eval(['save ' pf2 ' phlabel -append;']);
    
otherwise;
end;    




