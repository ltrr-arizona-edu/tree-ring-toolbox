function rwcomp
% rwcomp: compare a pair of measured ring-width series in measurement structure
% rwcomp;
% Last revised 10-23-02
%
% Compare a pair of measured ring-width series in a measurement structure
%
%*** INPUT
%
% No input arguments.
%
% User prompted to click on a .mat file that holds the measurement structure (see rwmeas)
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
% rwv option:  total ring width data exclusively taken from measured total ring width
% rwc option: total ring width data exclusively taken from computed (EW+LW) data (file XEL)
% rw option. total ring width data is taken from the partial-width data (EW+LW), or if that is not available, from the
%   measured total-width data (XT).  Thus the priority is for the partial-width measurements.  If data exist for both
%   total and partial width for a given year, priority is given to the partial-width data.

% check path1

% Get measurements

[file1,path1]=uigetfile('*.mat','input .mat storage file with the measurments');
pf1=[path1 file1];
eval(['load ' pf1 ';']);

kwh1=1;
while kwh1;
    kmen1=menu('Choose Mode of Comparison',...
        'Measured total ring width vs total ring width computed from EWW and LWW',...
        'Two measured total ring-width series',...
        'Two EWW series',...
        'Two LWW series',...
        'Two computed-total-width series');
    if kmen1==1; 
        kwh1=0;
        vbl2='XT';
        vbl1='XEL';
    elseif kmen1==2;  %  'Two measured total ring-width series',.
         kwh1=0;
         vbl2='XT';
         vbl1='XT';
     elseif kmen1==3; % two EWW series
         kwh1=0;
         vbl2='XEL';
         vbl1='XEL';
     elseif kmen1==4; % two LWW series
         kwh1=0;
         vbl2='XEL';
         vbl1='XEL';
     elseif kmen1==5; % Two computed-total-width series
         kwh1=0;
         vbl2='XEL';
         vbl1='XEL';
     else;
         uiwait(msgbox('First four modes only modes implemented so far','Click again!','modal'));
     end;
 end;
 
 
 kwh1=1;  
 
 while kwh1==1;
    kmen2=menu('Choose','Check Series','Quit');
    if kmen2==2; 
        kwh1=0;
    else;
        
        kwh2=1;
        while kwh2==1;
            close all;
            
            % Get partial width series
            if kmen1==1;
                nallow=[1 1];
            elseif kmen1==2;
                nallow=[2 2];
            elseif kmen1==3; % 2 EWW series
                nallow=[2 2];
            elseif kmen1==4; % 2 LWW series
                nallow=[2 2];
            elseif kmen1==5; % 2 computed-total-width series
                nallow=[2 2];
            else;
            end;
            eval(['nser1 = length(' vbl1 '.id);']);
            eval(['idset =' vbl1 '.id;']);
            nox1 = cellstr(repmat('-N',nser1,1));
            strmenu='Toggle N to Y to select series';       
            C1=idset; % cell array of names
            C2=nox1;
            Lin = logical(zeros(nser1,1));
            [C3,Lout] = menucell (C1,C2,Lin,nallow,strmenu); % C3 is cell array of selected series
            i1=find(Lout);
            
            % Check selected partial-width id against available total-width ids
            if kmen1==1;
                eval(['idset2 =' vbl2 '.id;']);
                Lwinner=strcmp(C3,idset2); 
            else;
            end;
            if kmen1==1 & sum(Lwinner)~=1;
                uiwait(msgbox(['No matching total-width series for ' char(C3)],...
                    'Pick again!','modal'));
            else
                if kmen1==1;
                    i2=find(Lwinner);
                else;
                    i2=i1(2);
                    i1=i1(1);
                end;
                
                if kmen1==1;
                    % Get the data , computed total from partials
                    eval(['x1 = ' vbl1 '.data{' int2str(i1) '}(:,4);']);
                    eval(['yrx1 = ' vbl1 '.data{' int2str(i1) '}(:,1);']);
                    % Get the measured total
                    eval(['x2 = ' vbl2 '.data{' int2str(i2) '}(:,2);']);
                    eval(['yrx2 = ' vbl2 '.data{' int2str(i2) '}(:,1);']);
                elseif kmen1==2;
                    % Get the data for first series
                    eval(['x1 = ' vbl1 '.data{' int2str(i1) '}(:,2);']);
                    eval(['yrx1 = ' vbl1 '.data{' int2str(i1) '}(:,1);']);
                    % Get the data for second series
                    eval(['x2 = ' vbl2 '.data{' int2str(i2) '}(:,2);']);
                    eval(['yrx2 = ' vbl2 '.data{' int2str(i2) '}(:,1);']);
                elseif kmen1==3; % two EWW
                    % Get the data for first series
                    eval(['x1 = ' vbl1 '.data{' int2str(i1) '}(:,2);']);
                    eval(['yrx1 = ' vbl1 '.data{' int2str(i1) '}(:,1);']);
                    % Get the data for second series
                    eval(['x2 = ' vbl2 '.data{' int2str(i2) '}(:,2);']);
                    eval(['yrx2 = ' vbl2 '.data{' int2str(i2) '}(:,1);']);
                elseif kmen1==4; % two LWW
                    % Get the data for first series
                    eval(['x1 = ' vbl1 '.data{' int2str(i1) '}(:,3);']);
                    eval(['yrx1 = ' vbl1 '.data{' int2str(i1) '}(:,1);']);
                    % Get the data for second series
                    eval(['x2 = ' vbl2 '.data{' int2str(i2) '}(:,3);']);
                    eval(['yrx2 = ' vbl2 '.data{' int2str(i2) '}(:,1);']);
                elseif kmen1==5;  % two computed-total-width series
                    % Get the data for first series
                    eval(['x1 = ' vbl1 '.data{' int2str(i1) '}(:,4);']);
                    eval(['yrx1 = ' vbl1 '.data{' int2str(i1) '}(:,1);']);
                    % Get the data for second series
                    eval(['x2 = ' vbl2 '.data{' int2str(i2) '}(:,4);']);
                    eval(['yrx2 = ' vbl2 '.data{' int2str(i2) '}(:,1);']);
                end;
                
                
                % Compute stats
                zmean =[nanmean(x1) nanmean(x2)]; % means of the two series
                zstd = [nanstd(x1) nanstd(x2)]; % std devs
                
                yrgo =max([yrx1(1) yrx2(1)]);
                yrsp = min([yrx1(end) yrx2(end)]);
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
                A = [a1(js) a2(js)] % the comp and meas 
                str1 = [num2str(B(:,1),'%6d') repmat('  ',5,1) num2str(B(:,2),'%5.2f')]
                str2 = [num2str(A,'%5.2f  ')];
                if kmen1==1;
                    str3=char('Year diff.    Comp   Meas',[str1 repmat('    ',5,1)  str2]);
                elseif kmen1==2;
                    str3=char(['Year diff.   ' C3{1}  ' ' C3{2}],  [str1 repmat('    ',5,1)  str2]);
                elseif kmen1==3 | kmen1==4 | kmen1==5;
                    str3=char(['Year diff.   ' C3{1}  ' ' C3{2}],  [str1 repmat('    ',5,1)  str2]);
                end;
                rat1 = max(s)/zstd(2); % error ratio = ratio of max absolute difference to std dev of measured ringwidths
                rat2= smed/zstd(2); % ratio of median abs difference to std dev of measured ringwidths
                str4=['Error Ratios: max = ' num2str(rat1,'%5.2f')  ';  median = ' num2str(rat2,'%5.3f')];
                             
                figure(1);
                if kmen1==1;
                    hp1=plot(yrx1,x1,'-o',yrx2,x2);
                elseif kmen1==2;
                    hp1=plot(yrx1,x1,'-o',yrx2,x2,'-^');
                elseif kmen1==3;
                    hp1=plot(yrx1,x1,'-o',yrx2,x2,'-^');
                elseif kmen1==4 | kmen1==5;
                    hp1=plot(yrx1,x1,'-o',yrx2,x2,'-^');
                end;
               
                set(hp1(2),'LineWidth',2);
                xlims=get(gca,'XLim');
                ylims=get(gca,'YLim');
                text(xlims(1),ylims(1),str3,'HorizontalAlignment','Left','VerticalAlignment','Bottom');
                if kmen1==1;
                    legend('Computed Total Width','Measured Total Width');
                elseif kmen1==2 | kmen1==3 | kmen1==4 | kmen1==5;
                    legend(C3{1},C3{2});
                end;
                ylabel('Width (mm)');
                xlabel('Year');
                set(gcf,'Position',[ 1          29        1024         672]);

                grid on;
                zoom xon;
                if kmen1==1;
                    title([char(C3) ': stand dev = ' num2str(zstd(1),'%5.3f') ';  ' str4]);
                elseif kmen1==2;
                    title(['Total Ring Width (measured):  Ave stand dev = ' num2str((zstd(1)+zstd(2))/2,'%5.3f') ';  ' str4]);
                elseif kmen1==3;
                    title(['Earlywood Width (measured):  Ave stand dev = ' num2str((zstd(1)+zstd(2))/2,'%5.3f') ';  ' str4]);
                    elseif kmen1==4;
                    title(['Latewood Width (measured):  Ave stand dev = ' num2str((zstd(1)+zstd(2))/2,'%5.3f') ';  ' str4]);
                elseif kmen1==5; 
                    title(['Computed total width:  Ave stand dev = ' num2str((zstd(1)+zstd(2))/2,'%5.3f') ';  ' str4]);
                end;
                
         
                
                kwh2=0;
            end;
            
            
        end; % while kwh2
      
    end; % while kwh1
end;

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




















    


