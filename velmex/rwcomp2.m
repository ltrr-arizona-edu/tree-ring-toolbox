function rwcomp2
% rwcomp2: compare pairs of measured ring-width series from two different rwmeas measurement files
% rwcomp2;
% Last revised 5-16-03
%
% compare pairs of measured ring-width series from two different rwmeas measurement files
%
%*** INPUT
%
% No input arguments.
%
% User prompted to click on a the two .mat files that hold the measurements
%
%*** OUTPUT
%
% No arguments. 
%
% Plots comparing series are shown on screen
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED 
%
% trailnan -- removes trailing NaNs
% uniqcell -- finds unique series names
% writerw -- writes a .rw, .eww or .lww file
%
%
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% Can compare a pair of series of the same type; four possible types; EWW, LWW, measured TW, computed TW 
%
% In automated mode, skips the graphics but gives summary ascii file
%
% rwv option:  total ring width data exclusively taken from measured total ring width
% rwc option: total ring width data exclusively taken from computed (EW+LW) data (file XEL)
% rw option. total ring width data is taken from the partial-width data (EW+LW), or if that is not available, from the
%   measured total-width data (XT).  Thus the priority is for the partial-width measurements.  If data exist for both
%   total and partial width for a given year, priority is given to the partial-width data.

% check path1


%************  AUTOMATED?

batchmode = questdlg('Batch mode -- automatic comparison of pairs of series; no plots?');


%************  GET THE TWO MEASUREMENT FILES

[file3,path3]=uigetfile('*.mat','First mat measurement file');
pf3=[path3 file3];
eval(['load ' pf3 ' XT XEL;']);
XT1=XT;
XEL1=XEL;
[file3a,path3a]=uigetfile('*.mat',['Second .mat measurement file (first file is ' file3 ' )']);
pf3a=[path3a file3a];
eval(['load ' pf3a ' XT XEL;']);
XT2=XT;
XEL2=XEL;




%************  TYPE OF MEASUREMENT TO BE COMPARED

kwh1=1;
while kwh1;
    kmen1=menu('Choose type of measurements to be compared',...
        'EWW',...
        'LWW',...
        'Computed total width',...
        'Measured total width',...
        'Abort');
    if kmen1==1; % earlywood
        kwh1=0;
        vbl2='XEL2';
        vbl1='XEL1';
    elseif kmen1==2;  % latewood
        kwh1=0;
        vbl2='XEL2';
        vbl1='XEL1';
    elseif kmen1==3; % computed total
        kwh1=0;
        vbl2='XEL2';
        vbl1='XEL1';
    elseif kmen1==4; % measured total
        kwh1=0;
        vbl2='XT2';
        vbl1='XT1';
    else;
        clc;
        disp('Mission Aborted, as your requested');
        return;
        %uiwait(msgbox('First four modes only modes implemented so far','Click again!','modal'));
    end;
end;
XX=vbl1;
YY=vbl2;
if kmen1<4;
    ringtype='partial';
else;
    ringtype='total';
end;


%***************************  CHECK THAT BOTH FILES HAVE AT LEAST ONE SERIES OF DESIRED TYPE

nlenx=eval(['length(' XX '.id);']);
nleny=eval(['length(' YY '.id);']);
if nlenx==1;
    eval(['Lcheck= isempty(' XX '.id{1});']);
    if Lcheck;
        error([file3 ' has no ' ringtype '-width measurements']);
    end;
    eval(['Lcheck= isempty(' YY '.id{1});']);
    if Lcheck;
        error([file3a ' has no ' ringtype '-width measurements']);
    end;
end;
msgtemp={[file3 ' has ' num2str(nlenx) ' '  ringtype '-width measurements'],...
        [file3a ' has ' num2str(nleny) ' '  ringtype '-width measurements']};
uiwait(msgbox(msgtemp,'Message','modal'));



%****************   SPECIFY HOW NAMES SHOULD MATCH

if strcmp(batchmode,'Yes');
    
    aliascode=('Map site code from one file to different code in another');
    
    % Find all unique 3-char opening codes for series in first file
    codes1=char(eval([XX '.id']));
    codes1=unique(cellstr(codes1(:,1:3)));
    ncodes1=length(codes1);
    
    % Find all unique 3-char opening codes for series in second file
    codes2=char(eval([YY '.id']));
    codes2=unique(cellstr(codes2(:,1:3)));
    ncodes2=length(codes2);
    
    
    
    
    
    matchmode=questdlg('Match all series with identical core ids?');
    if ~strcmp(matchmode,'Yes');
        prompt={['Match series from ' file3 ' with 3-letter code: '],['with series from ' file3a ' with 3-letter code: ']};
        def={'',''};
        dlgTitle='ID-Matching';
        lineNo=1;
        answer=inputdlg(prompt,dlgTitle,lineNo,def);
        idcode1=upper(answer{1});
        idcode2=upper(answer{2});
    end; % if ~strcmp(matchmode,'Yes');
else; % if strcmp(batchmode,'Yes');
end; % if strcmp(batchmode,'Yes');


if ~strcmp(batchmode,'Yes'); %  This section only for manual comparison of a pair of series
    
    
    
    kwh1=1;  
    
    while kwh1==1;
        kmen2=menu('Choose','Compare a pair of series','Quit');
        if kmen2==2; 
            kwh1=0;
        else;
            
            kwh2=1;
            while kwh2==1;
                close all;
                
                nallow=[1 1]; % max of one series to be selected from each of the two files
                
                % Get a series from first file
                nser1=nlenx; % number of series in file
                eval(['idset =' XX '.id;']);
                nox1 = cellstr(repmat('-N',nser1,1));
                strmenu=['Toggle N to Y to select a series from ' file3];       
                C1=idset; % cell array of names
                C2=nox1;
                Lin = logical(zeros(nser1,1));
                [C3,Lout] = menucell (C1,C2,Lin,nallow,strmenu); % C3 is cell array of selected series
                i1=find(Lout);
                
                % Get a series from second file  ( note the 2 in nser2, and y as last char in names)
                nser2=nleny; % number of series in file
                eval(['idsety =' YY '.id;']);
                nox1y= cellstr(repmat('-N',nser2,1));
                strmenuy=['Toggle N to Y to select a series from ' file3a ': MATCHING ' C3{1}];       
                C1y=idsety; % cell array of names
                C2y=nox1y;
                Liny = logical(zeros(nser2,1));
                [C3y,Louty] = menucell (C1y,C2y,Liny,nallow,strmenuy); % C3y is cell array of selected series
                i1y=find(Louty);
                i2=i1y; % rename
                
                
                %***************** GET THE DATA FOR COMPARISON
                if kmen1==1; % eww
                    % Get first
                    eval(['x1 = ' XX '.data{' int2str(i1) '}(:,2);']);
                    eval(['yrx1 = ' XX '.data{' int2str(i1) '}(:,1);']);
                    % Get second
                    eval(['x2 = ' YY '.data{' int2str(i2) '}(:,2);']);
                    eval(['yrx2 = ' YY '.data{' int2str(i2) '}(:,1);']);
                elseif kmen1==2; % LWW
                    % Get the data for first series
                    eval(['x1 = ' XX '.data{' int2str(i1) '}(:,3);']);
                    eval(['yrx1 = ' XX '.data{' int2str(i1) '}(:,1);']);
                    % Get the data for second series
                    eval(['x2 = ' YY '.data{' int2str(i2) '}(:,3);']);
                    eval(['yrx2 = ' YY '.data{' int2str(i2) '}(:,1);']);
                elseif kmen1==3; % total width computed from partial
                    % Get the data for first series
                    eval(['x1 = ' XX '.data{' int2str(i1) '}(:,4);']);
                    eval(['yrx1 = ' XX '.data{' int2str(i1) '}(:,1);']);
                    % Get the data for second series
                    eval(['x2 = ' YY '.data{' int2str(i2) '}(:,4);']);
                    eval(['yrx2 = ' YY '.data{' int2str(i2) '}(:,1);']);
                elseif kmen1==4; % measured total width 
                    % Get the data for first series
                    eval(['x1 = ' XX '.data{' int2str(i1) '}(:,2);']);
                    eval(['yrx1 = ' XX '.data{' int2str(i1) '}(:,1);']);
                    % Get the data for second series
                    eval(['x2 = ' YY '.data{' int2str(i2) '}(:,2);']);
                    eval(['yrx2 = ' YY '.data{' int2str(i2) '}(:,1);']);
                end;
                
                % Compute stats
                zmean =[nanmean(x1) nanmean(x2)]; % means of the two series
                zstd = [nanstd(x1) nanstd(x2)]; % std devs
                
                yrgo =max([yrx1(1) yrx2(1)]);
                yrsp = min([yrx1(end) yrx2(end)]);
                
                if (yrsp-yrgo)<3;
                    error(['Fewer than 3 years overlap in ' C3{1} ' from ' file3 ' and ' C3y{1} ' from file3a']);
                end;
                
                L1 = yrx1>=yrgo & yrx1<=yrsp;
                L2 = yrx2>=yrgo & yrx2<=yrsp;
                a1 = x1(L1);
                yra1=yrx1(L1);
                a2=x2(L2);
                yra2=yrx2(L2);
                d = [a2 - a1]; % difference measured minus computed, or second minus first
                [s js] = sort(abs(d)); % sort low to high max absolute differences
                smed = nanmedian(s); % median absolute error
                % truncate NaN rows if any
                if any(isnan(s));
                    Lnan =isnan(s);
                    ibad =find(Lnan);
                    s(ibad(1):end)=[];
                    js(ibad(1):end)=[];
                end;
                % want highest five
                s=flipud(s);
                js=flipud(js);
                s=s(1:5);
                js=js(1:5);
                B = [yra2(js) s];
                A = [a1(js) a2(js)] ;  % the comp and meas 
                str1 = [num2str(B(:,1),'%6d') repmat('  ',5,1) num2str(B(:,2),'%5.2f')];
                str2 = [num2str(A,'%5.2f  ')];
                str3=char(['Year diff.   ' C3{1}  ' ' C3y{1}],  [str1 repmat('    ',5,1)  str2]);
                rat1 = max(s)/zstd(2); % error ratio = ratio of max absolute difference to std dev of measured ringwidths
                rat2= smed/zstd(2); % ratio of median abs difference to std dev of measured ringwidths
                str4=['Error Ratios: max = ' num2str(rat1,'%5.2f')  ';  median = ' num2str(rat2,'%5.3f')];
                
                
                %------------ GRAPHICS
                
                figure(1);
                hp1=plot(yrx1,x1,'-o',yrx2,x2,'-^');
                
                set(hp1(2),'LineWidth',2);
                xlims=get(gca,'XLim');
                ylims=get(gca,'YLim');
                text(xlims(1),ylims(1),str3,'HorizontalAlignment','Left','VerticalAlignment','Bottom');
                
                legend([C3{1} ', from ' file3],[C3y{1} ', from ' file3a]);
                
                
                ylabel('Width (mm)');
                xlabel('Year');
                set(gcf,'Position',[ 1          29        1024         672]);
                
                grid on;
                zoom xon;
                
                % Build title
                if kmen1==3;
                    txta = 'Computed ';
                else;
                    txta='Measured ';
                end;
                if kmen1==1;
                    txtb='EWW';
                elseif kmen1==2;
                    txtb = 'LWW';
                else;
                    txtb='Total Width';
                end;
                txtc=[txta txtb];
                title([txtc ':  Ave stand dev = ' num2str((zstd(1)+zstd(2))/2,'%5.3f') ';  ' str4]);
                
                kwh2=0;
            end; % if kmen2=2
        end; % while kwh2
    end% while kwh1
    
    
    
    
    
    
    
    
    
    
else; %  if strcmp(batchmode,'Yes'); %  This section only for manual comparison of a pair of series; HERE CODE FOR AUTOMATED 
end; % if strcmp(batchmode,'Yes'); %  This section only for manual comparison of a pair of series


%  %--------------------  OLD CODE  
%  
%             % Check selected partial-width id against available total-width ids
%             if kmen1==1;
%                 eval(['idset2 =' vbl2 '.id;']);
%                 Lwinner=strcmp(C3,idset2); 
%             else;
%             end;
%             if kmen1==1 & sum(Lwinner)~=1;
%                 uiwait(msgbox(['No matching total-width series for ' char(C3)],...
%                     'Pick again!','modal'));
%             else
%                 if kmen1==1;
%                     i2=find(Lwinner);
%                 else;
%                     i2=i1(2);
%                     i1=i1(1);
%                 end;
%                 
%                 if kmen1==1;
%                     % Get the data , computed total from partials
%                     eval(['x1 = ' vbl1 '.data{' int2str(i1) '}(:,4);']);
%                     eval(['yrx1 = ' vbl1 '.data{' int2str(i1) '}(:,1);']);
%                     % Get the measured total
%                     eval(['x2 = ' vbl2 '.data{' int2str(i2) '}(:,2);']);
%                     eval(['yrx2 = ' vbl2 '.data{' int2str(i2) '}(:,1);']);
%                 elseif kmen1==2;
%                     % Get the data for first series
%                     eval(['x1 = ' vbl1 '.data{' int2str(i1) '}(:,2);']);
%                     eval(['yrx1 = ' vbl1 '.data{' int2str(i1) '}(:,1);']);
%                     % Get the data for second series
%                     eval(['x2 = ' vbl2 '.data{' int2str(i2) '}(:,2);']);
%                     eval(['yrx2 = ' vbl2 '.data{' int2str(i2) '}(:,1);']);
%                 elseif kmen1==3; % two EWW
%                     % Get the data for first series
%                     eval(['x1 = ' vbl1 '.data{' int2str(i1) '}(:,2);']);
%                     eval(['yrx1 = ' vbl1 '.data{' int2str(i1) '}(:,1);']);
%                     % Get the data for second series
%                     eval(['x2 = ' vbl2 '.data{' int2str(i2) '}(:,2);']);
%                     eval(['yrx2 = ' vbl2 '.data{' int2str(i2) '}(:,1);']);
%                 elseif kmen1==4; % two LWW
%                     % Get the data for first series
%                     eval(['x1 = ' vbl1 '.data{' int2str(i1) '}(:,3);']);
%                     eval(['yrx1 = ' vbl1 '.data{' int2str(i1) '}(:,1);']);
%                     % Get the data for second series
%                     eval(['x2 = ' vbl2 '.data{' int2str(i2) '}(:,3);']);
%                     eval(['yrx2 = ' vbl2 '.data{' int2str(i2) '}(:,1);']);
%                 elseif kmen1==5;  % two computed-total-width series
%                     % Get the data for first series
%                     eval(['x1 = ' vbl1 '.data{' int2str(i1) '}(:,4);']);
%                     eval(['yrx1 = ' vbl1 '.data{' int2str(i1) '}(:,1);']);
%                     % Get the data for second series
%                     eval(['x2 = ' vbl2 '.data{' int2str(i2) '}(:,4);']);
%                     eval(['yrx2 = ' vbl2 '.data{' int2str(i2) '}(:,1);']);
%                 end;
%                 
%                 
%                 
%                 
%                 figure(1);
%                 if kmen1==1;
%                     hp1=plot(yrx1,x1,'-o',yrx2,x2);
%                 elseif kmen1==2;
%                     hp1=plot(yrx1,x1,'-o',yrx2,x2,'-^');
%                 elseif kmen1==3;
%                     hp1=plot(yrx1,x1,'-o',yrx2,x2,'-^');
%                 elseif kmen1==4 | kmen1==5;
%                     hp1=plot(yrx1,x1,'-o',yrx2,x2,'-^');
%                 end;
%                
%                 set(hp1(2),'LineWidth',2);
%                 xlims=get(gca,'XLim');
%                 ylims=get(gca,'YLim');
%                 text(xlims(1),ylims(1),str3,'HorizontalAlignment','Left','VerticalAlignment','Bottom');
%                 if kmen1==1;
%                     legend('Computed Total Width','Measured Total Width');
%                 elseif kmen1==2 | kmen1==3 | kmen1==4 | kmen1==5;
%                     legend(C3{1},C3{2});
%                 end;
%                 ylabel('Width (mm)');
%                 xlabel('Year');
%                 set(gcf,'Position',[ 1          29        1024         672]);
% 
%                 grid on;
%                 zoom xon;
%                 if kmen1==1;
%                     title([char(C3) ': stand dev = ' num2str(zstd(1),'%5.3f') ';  ' str4]);
%                 elseif kmen1==2;
%                     title(['Total Ring Width (measured):  Ave stand dev = ' num2str((zstd(1)+zstd(2))/2,'%5.3f') ';  ' str4]);
%                 elseif kmen1==3;
%                     title(['Earlywood Width (measured):  Ave stand dev = ' num2str((zstd(1)+zstd(2))/2,'%5.3f') ';  ' str4]);
%                     elseif kmen1==4;
%                     title(['Latewood Width (measured):  Ave stand dev = ' num2str((zstd(1)+zstd(2))/2,'%5.3f') ';  ' str4]);
%                 elseif kmen1==5; 
%                     title(['Computed total width:  Ave stand dev = ' num2str((zstd(1)+zstd(2))/2,'%5.3f') ';  ' str4]);
%                 end;
%                 
%          
%                 
%                 kwh2=0;
%             end;
%             
%          
% end;



quest1=questdlg('Close Windows');
if strcmp(quest1,'Yes');
    close all;
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























