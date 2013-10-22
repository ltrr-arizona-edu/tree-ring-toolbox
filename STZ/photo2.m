function S=photo2



%*** NOTES
%
% cd to directory with photo files before running


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
T=S;

% While loop over files
kwh1=1; 
while kwh1;
    T1=sort(T(1:nfile1));
    T=T1;
    T{nfile1+1}='Quit';
    for k=1:3;
        if ishandle(k);
            close(k);
        end;
    end;
        
    kmen1 = menu('Choose ',T);
    if kmen1==nfile1+1;
        kwh1=0;
    else;
        kwh1=1;
        f = T{kmen1};
        a=imread(f);
        aorig=a;
        figure(1);
        imshow(aorig);
        title([f ': original']);
        figure(2); 
        imshow(a);
        porig=get(gcf,'Position'); % figure position initial
        title([f ': current version']);
        
        % While loop to operate on files
        kwh2=1;
        while kwh2;
            kmen2=menu('Choose',...
                'Delete file from disk',...
                'Rotate',...
                'Crop',...
                'Resize',...
                'RGB2gray',...
                'brighten',...
                'Adjust RGB',...
                'Label',...
                'Replace original file with version in Figure 2',...
                'Save special',...
                'Add to current subplot',...
                'Return to file menu');
            switch kmen2;
            case 1; % delete file
                kquest = questdlg(['Are you sure you want to delete ' f ' from disk?']);
                if strcmp(kquest,'Yes');
                    
                    eval(['delete ' f ';']);
                    T{kmen1}=['Deleted: '  f];
                else;
                end;
            case 2; % rotate;
                kfirst=1;
                atemp=a;
                kwh3=1; 
                if kfirst==1;
                    krot=1;
                else;
                    krot=NaN;
                end;
                while kwh3;
                    if kfirst~=1;
                        krot=menu('Choose',...
                            'Rotate again',...
                            'Revert to image before rotation',...
                            'Accept rotation and go up a menu');
                    else;
                        kfirst=0;
                        krot=1;
                    end;
                    
                    switch krot;
                    case 3; % Accept rotation
                        kwh3=0;
                        a=atemp;
                        close(3);
                        figure(2);
                        imshow(a);
                        break;
                    case 2; % revert to image before rotation 
                        kwh3=1;
                        atemp=a;
                        figure(3);
                        imshow(atemp);
                        title([f ' reverted to before rotation image']);
                    case 1; % rotate first time or rotate again again
                        kwh3=1;
                                      
                        % Positive is CCW
                        prompt={'Enter degrees to rotate:'};
                        def={'90'};
                        dlgTitle='Rotation (+ is CCW)';
                        lineNo=1;
                        answer=inputdlg(prompt,dlgTitle,lineNo,def);
                        degrot=str2num(answer{1});
                        atemp=imrotate(atemp,degrot);
                        figure(3);
                        imshow(atemp);
                        title([f ' rotated ']);
                    end; % switch krot
                end; % while kwh3
                
                
            case 3; % Crop
                kfirst=1;
                atemp=a;
                kwh3=1; 
                if kfirst==1;
                    kcrop=1;
                else;
                    kcrop=NaN;
                end;
                while kwh3;
                    if kfirst~=1;
                        kcrop=menu('Choose',...
                            'Crop more',...
                            'Revert to image before cropping',...
                            'Accept cropping and go up a menu');
                    else;
                        kfirst=0;
                        kcrop=1;
                    end;
                    
                    switch kcrop;
                    case 3; % Accept cropping
                        kwh3=0;
                        a=atemp;
                        close(3);
                        figure(2);
                        imshow(a);
                        break;
                    case 2; % revert to image before cropping 
                        kwh3=1;
                        atemp=a;
                        figure(3);
                        imshow(atemp);
                        title([f ' reverted to before-cropping image']);
                    case 1; % crop first time or crop again
                        kwh3=1;
                        atemp=imcrop(atemp);
                        figure(3);
                        imshow(atemp);
                        title([f ' cropped ']);
                    end; % switch kcrop
                end; % while kwh3
                
            case 4; % resize image matrix
                kfirst=1;
                atemp=a;
                kwh3=1; 
                if kfirst==1;
                    kresize=1;
                else;
                    kresize=NaN;
                end;
                while kwh3;
                    if kfirst~=1;
                        kresize=menu('Choose',...
                            'Resize again',...
                            'Revert to image before resizing',...
                            'Accept resizing and go up a menu');
                    else;
                        
                        kresize=1;
                    end;
                    
                    switch kresize;
                    case 1; % resize first time or resize again
                        if kfirst==1;
                            fcurrent='1';
                            methcurrent='bilinear';
                        end;
                        
                        prompt={'Enter resizing factor','Enter Method (nearest,bilinear,bicubic'};
                        
                        def={fcurrent,methcurrent};
                        dlgTitle='Resizing factor (<1 decreases resolution)';
                        lineNo=1;
                        answer=inputdlg(prompt,dlgTitle,lineNo,def);
                        xresize=str2num(answer{1});
                        xmeth =answer{2};
                        fcurrent=num2str(xresize);
                        methcurrent=xmeth;
                        atemp=imresize(atemp,xresize,xmeth);
                        figure(3);
                        imshow(atemp);
                        title([f ' resized at x ' num2str(xresize)]);
                        set(gcf,'Position',porig);
                        kfirst=0;
                        kwh3=1;
                        
                    case 2; % revert to image before cropping 
                        atemp=a;
                        figure(3);
                        imshow(atemp);
                        title([f ' reverted to before-resizing image']);
                        kfirst=1;
                        kwh3=1;
                    case 3; % Accept resizing
                        kwh3=0;
                        a=atemp;
                        close(3);
                        figure(2);
                        imshow(a);
                        set(gcf,'Position',porig);
                        break;
                        
                    end; % switch kresize
                end; % while kwh3
                
                
                
            case 5; % convert to grayscale
                kfirst=1;
                atemp=a; % store current version of image
                kwh3=1; % while control
                while kwh3;
                    if kfirst~=1;
                        kmen3= menu('Choose one',...
                            'Toggle rgb - grayscale',...
                            'Accept current version -- gray or rgb');
                        switch kmen3;
                        case 1; % revert to rgb
                            kwh3=1;
                            if ~isgray(atemp);
                                atemp=rgb2gray(atemp);
                            else;
                                atemp=a;
                            end;
                            figure(3);
                            imshow(atemp);
                        case 2; % accept revision
                            kwh3=0; 
                            a=atemp;
                            close(3);
                            figure(2);
                            imshow(a);
                            title([f ' -- after grayscale conversion']);
                        end; % switch kmen3
                    else; % first time -- toggle
                        kwh3=1;
                        kfirst=0;
                        if isgray(atemp);
                            error('image already grayscale');
                        else;
                            atemp=rgb2gray(atemp);
                            figure(3);
                            imshow(atemp);
                            title([f ': grayscaled']);
                        end; % if kfirst
                    end; % while kwh3
                end; % wwhile kwh3
                
                
            case 6; % brighten
                kfirst=1;
                atemp=a;
                kwh3=1; 
                definit={'1'};
                defprev=definit;
                if kfirst==1;
                    kbright=1;
                else;
                    kbright=NaN;
                end;
                while kwh3;
                    if kfirst~=1;
                        kbright=menu('Choose',...
                            'Brighten again',...
                            'Revert to image before brightening',...
                            'Accept brightening and go up a menu');
                    else;
                    end;
                    switch kbright;
                    case 1; % brighten first time or brighten again
                        kwh3=1;
                        % gamma>1 weights towards darker
                        prompt={'Enter brighten parameter:'};
                        if kfirst==1;
                            def=definit;
                        else;
                            def=defprev;
                        end;
                        dlgTitle='Brightening (gamma>1 will darken)';
                        lineNo=1;
                        answer=inputdlg(prompt,dlgTitle,lineNo,def);
                        defprev={answer{1}};
                        gamma=str2num(answer{1});
                        atemp = imadjust(a,[0 0 0; 1 1 1],[0 0 0; 1 1 1],gamma);
                        figure(3);
                        imshow(atemp);
                        title([f ' brightened with gamma = ' num2str(gamma)]);
                        kfirst=0;
                    case 2; % revert to image before brightening 
                        kwh3=1;
                        atemp=a;
                        figure(3);
                        imshow(atemp);
                        title([f ' reverted to before brightening image']);
                    case 3; % Accept brightening
                        kwh3=0;
                        a=atemp;
                        close(3);
                        figure(2);
                        imshow(a);
                        break;
                        
                        
                    end; % switch kbright
                end; % while kwh3
                
            case 7; % Adjust rgb
                definit={'[0 0 0]','[1 1 1]','[0 0 0]','[1 1 .85]','[1 1 1]'};
                defprev=definit;
                kfirst=1;
                atemp=a;
                kwh3=1; 
                if kfirst==1;
                    kadjust=1;
                else;
                    kadjust=NaN;
                end;
                while kwh3;
                    if kfirst~=1;
                        kadjust=menu('Choose',...
                            'Adjust again',...
                            'Revert to image before adjusting',...
                            'Accept adjusting and go up a menu');
                    else;
                    end;
                    switch kadjust;
                    case 1; % adjust first time or again
                        kwh3=1;
                        prompt={'Enter Lowin :', 'Enter High in','Enter Lo-out','Enter Hi-out','Enter gamma'};
                        if kfirst==1;
                            def=definit;
                        else;
                            def=defprev;
                        end;
                        dlgTitle='RGB Adjustment)';
                        lineNo=1;
                        answer=inputdlg(prompt,dlgTitle,lineNo,def);
                        lowin=str2num(answer{1});
                        highin=str2num(answer{2});
                        lowout=str2num(answer{3});
                        highout=str2num(answer{4});
                        gamma=str2num(answer{5});
                        defprev={answer{1},answer{2},answer{3},answer{4},answer{5}};
                                               
                        atemp=imadjust(a,[lowin; highin],[lowout; highout],gamma);
                        figure(3);
                        imshow(atemp);
                        title([f ' adjusted']);
                        kfirst=0;
                        
                    case 2; % revert to image before adjusting 
                        kwh3=1;
                        atemp=a;
                        figure(3);
                        imshow(atemp);
                        title([f ' reverted to before-adjusted image']);
                    case 3; % Accept adjusting
                        kwh3=0;
                        a=atemp;
                        close(3);
                        figure(2);
                        imshow(a);
                        break;
                    end; % switch kadjust
                end; % while kwh3
                
            case 8; % Label
                
                kfirst=1;
                atemp=a;
                kwh3=1; 
                if kfirst==1;
                    klabel=1;
                else;
                    klabel=NaN;
                end;
                while kwh3;
                    if kfirst~=1;
                        klabel=menu('Choose',...
                            'Label again',...
                            'Revert to image before labeling',...
                            'Accept labeling and go up a menu');
                    else;
                        kfirst=0;
                        klabel=1;
                    end;
                    switch klabel;
                    case 3; % Accept labeling
                        kwh3=0;
                        
                        break;
                    case 2; % revert to image before labeling 
                        kwh3=1;
                        atemp=a;
                        figure(2);
                        imshow(atemp);
                        title([f ' reverted to before-labeling image']);
                    case 1; % label first time or label again
                        kwh3=1;
                        figure(2);
                        prompt={'Enter Font Size:','Enter Label:'};
                        fdef = f;
                        def={'20',fdef};
                        dlgTitle='Label Data';
                        lineNo=1;
                        answer=inputdlg(prompt,dlgTitle,lineNo,def);
                        fs1 = str2num(answer{1});
                        txt1 = char(answer{2});
                        gtext(txt1,'FontSize',fs1);
                        
                    end; % switch klabel
                end; % while kwh3
                
            case 9; % replace original file
                % Prompt for filename
                prompt={'Enter name:'};
                def={f};
                dlgTitle='New name to imwrite file to';
                lineNo=1;
                answer=inputdlg(prompt,dlgTitle,lineNo,def);
                fout=answer{1};
                
                % Prompt for quality
                kqual=menu('Desired quality value of .jpg file',...
                    '25',...
                    '50',...
                    '75');
                switch kqual;
                case 1;
                    kq = 25;
                case 2; 
                    kq = 50;
                case 3;
                    kq=75;
                end; %switc kqual
                
                % Write, checking that file with same name as writefile does not exist
                L=exist(fout,'file')==2;
                if L; % file of same name already exists
                    kquest=questdlg(['Sure you want to overwrite exising ' fout '?']);
                    if strcmp(kquest,'Yes');
                        imwrite(a,fout,'JPG','Quality',kq);
                        T{kmen1}=['Overwritten: ' fout];
                    else; 
                    end;
                else; % file of same name does note yet exist
                    imwrite(a,fout,'JPG','Quality',kq);
                    T{kmen1}=['Renamed: ' fout];
                    % Delete old file
                    kquest = questdlg(['Are you sure you want to delete ' f ' from disk?']);
                    if strcmp(kquest,'Yes');
                        
                        eval(['delete ' f ';']);
                        
                    else;
                    end;
                    
                end;
               
            case 10; % save special
                
                % Prompt for filename
                prompt={'Enter name:'};
                def={['s-' f]};
                dlgTitle='Name to save special version as';
                lineNo=1;
                answer=inputdlg(prompt,dlgTitle,lineNo,def);
                fout=answer{1};
                
                 % Prompt for quality
                kqual=menu('Desired quality value of .jpg file',...
                    '25',...
                    '50',...
                    '75');
                switch kqual;
                case 1;
                    kq = 25;
                case 2; 
                    kq = 50;
                case 3;
                    kq=75;
                end;
                
                % Write
                imwrite(a,fout,'JPG','Quality',kq);
                
            case 11; % add to current subplot
                figure(2)
                jsub=jsub+1;
                
                % Initialize subplot
                if jsub==1;
                    
                    prompt={'Enter number of rows:','Enter number of columns:','Enter number of cells to fill:'};
                    def={'3','2','6'};
                    dlgTitle='Dimensions of subplot';
                    lineNo=1;
                    answer=inputdlg(prompt,dlgTitle,lineNo,def);
                    nrow=str2num(answer{1});
                    ncol=str2num(answer{2});
                    ncell=str2num(answer{3});
                    nsub=nrow*ncol;
                    shrink=(  round(min([ncol nrow]))  )/ nsub;
                    
                    phcolor=cell(ncell,1);
                    phname=cell(ncell,1);
                    phlabel=cell(ncell,1);
                    paxis=repmat(NaN,ncell,4);
                    scalef=repmat(NaN,ncell,1);
                    xylabel=repmat(NaN,ncell,2);
                    
                    
                    
                end;
                atemp=imresize(a,shrink,'bilinear');
                                   
                
                
                if jsub>ncell;
                    jsub=0;
                    nwind=nwind+1;
                end;
                figure(nwind);
                set(gcf,'Position',[232   258   315   420],'Color',[1 1 1]);
                set(gcf,'PaperPosition',[.25 1 6 8]);
                subplot(nrow,ncol,jsub);
                subimage(atemp);
                set(gca,'XTick',[],'YTick',[]);
                
                % Scale image, depending whether in left or right column 
                p1=get(gca,'Position');
                
                % Initial scaling expands subplots to reduce white space
                scfact=1.3; % 
                [p1,scfact]=subfcn2(p1,scfact);
                set(gca,'Position',p1);
                                        
                                
                % While loop for rescaling subplot
                kwhscale=1;
                while kwhscale;
                    kquest=questdlg('Rescale subplot?');
                    switch kquest;
                    case 'Yes';
                        
                        prompt={'Enter scaling factor:'};
                        def={'1.0'};
                        dlgTitle='Optional re-scaling of subplot';
                        lineNo=1;
                        answer=inputdlg(prompt,dlgTitle,lineNo,def);
                        scfact=str2num(answer{1});
                        [p1,scfact]=subfcn2(p1,scfact);
                        kwhscale=1;
                        set(gca,'Position',p1);
                    otherwise;
                        kwhscale=0;
                    end; % swithch kquest
                    
                    
                end; % while kwhscale
                
                
                % LABELING
                
                kfirst=1;
                atemp=a;
                kwh3=1; 
                if kfirst==1;
                    klabel=1;
                else;
                    klabel=NaN;
                end;
                while kwh3;
                    if kfirst~=1;
                        klabel=menu('Choose',...
                            'Label again',...
                            'Revert to image before labeling',...
                            'Accept labeling and go up a menu');
                    else;
                        %kfirst=0;
                        klabel=1;
                    end;
                    switch klabel;
                    case 3; % Accept labeling
                        kwh3=0;
                        
                        phname{jsub}=f;
                        scalef(jsub)=scfact;
                        paxis(jsub,:) = get(gca,'Position');
                        phlabel{jsub}=get(tt,'String');
                        xytemp=get(tt,'Position');
                        xylabel(jsub,:)=xytemp(1:2);
                        if isgray(atemp);
                            phcolor{jsub}='Gray';
                        else;
                            phcolor{jsub}='RGB';
                        end;
                        
                        % Produce .fig file
                        if jsub==ncell;
                            pfig=get(gcf,'Position');
                            pappos = get(gcf,'PaperPosition');
                            fslab=fs1;
                            pagenm = ['page' int2str(nwind-3) '.fig'];
                            pf2=uiputfile(pagenm,'Outfile (.fig) for figure of subplots');
                            eval(['hgsave ' pf2 ';']);
                            matty=strtok(pf2,'.');
                            
                            set1=[' pagenm nrow ncol ncell pappos  pfig fslab shrink '];
                            set2=[' phname scalef phcolor paxis phlabel xylabel '];  
                            eval(['save ' matty set1 set2 ';'])
                            
                            jsub=0;
                            
                        end;
                        
                        break;
                    case 2; % revert to image before labeling 
                        kwh3=1;
                        %atemp=a;
                        figure(nwind);
                        imshow(atemp);
                        %title([f ' reverted to before-labeling image']);
                    case 1; % label first time or label again
                        kwh3=1;
                        figure(nwind);
                        x=get(gca,'XLim');
                        y=get(gca,'YLim');
                        fdef = strtok(f,'.');
                        if kfirst~=1;
                            prompt={'Enter Font Size:','Enter Label:'};
                            def={'10',fdef};
                            dlgTitle='Label Data';
                            lineNo=1;
                            answer=inputdlg(prompt,dlgTitle,lineNo,def);
                            fs1 = str2num(answer{1});
                            txt1 = char(answer{2});
                            set(tt,'String',txt1,'FontSize',fs1);
                            %text(x(1),y(2),txt1,'VerticalAlignment','Top')
                        else; % kfirst==1
                            
                            fs1=10;
                            tt=text(x(1),y(2),fdef,'VerticalAlignment','Top')
                            kfirst=0;
                        end;
                        
                    end; % switch klabel
                end; % while kwh3
                
                
            case 12; % Return fo file menu
                kwh2=0;
                
            end; % switch kmen2
            
            
            
        end; % while kwh2
    end;
end;





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








