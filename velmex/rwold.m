function rwmeas
% rwmeas:  ring-width measurement
% rwmeas;
% Last revised 9-24-01
%
% Measure ring widths with Velmex system.  Options for earlywood/latewood width or total width
%
%*** INPUT
%
% No args.  User prompted for files and settings
%
%
%*** OUTPUT 
%
% No args. User prompted for filenames.  File created or modified is .mat storage of 3-d matrix
% with measured series.  
%
%
%*** NOTES
%
% Can append data to .mat file, or create new file
% Use other functions ?? for converting data to .rw or .rwl
% Data to be stored in cells within structure:
%   XT, XEL  -- structures 
%      Fields:
%       .units -- string measurement units (e.g., 'mm x 100')
%       .summary -- char array of 1)id, 2)span measured total rw, 3-5) spans of EW,LW, sum of EW+LW
%       NEXT FIELDS CONTAIN CELLS, WITH ONE ELEMENT PER SERIES
%       .data -- cell with tsm for each measured series, year as col 1. (2 col for total, 4 col for EL)
%           For example, XT.data{1} has a tsm of total width for series
%       .id --- cell with id of each series
%       .span -- first, last year of series 
%       .who -- measurer -- also a cell variable
%       .when -- date measurement completed --- also in cells
%
% For example, XEL.data{3} is a 4-column tsm of year , EWW, LWW, EWW+LWW for series 3. 
% XEL.span{3} is a  3x2 matrix with first and last years (span) of EWW, LWW, EWW+LWW 
% XEL.who{3} holds a string variable of who measured --- e.g.,   'DMM'
% XEL.when{3} holds date and time last measurements on the series made -- e.g., 08-05-01
% XEL.data{} has 4 columns: year, EWW, LWW and total of EWW+LWW.  This last column of total ring width contrasts with
%   the contents of XT, which is directly measured total ring width.
%
% Series are sorted and stored in alphabetical order by id

% % Hard Code
% units='mm x 100';
% person='DMM';
clear;
close all;
clc;
set1=' XT XEL vlist';


%---- MEASURE OR EDIT?

kmen=menu('What doe you want to do?','Measure','Edit');
if kmen==1;
    runmode='Measure';
else;
    runmode='Edit';
end;
clear kmen;


%--- CREATE NEW SITE OR OPERATE ON EXISTING SITE

switch runmode;
case 'Measure';
    
case 'Edit';
    filemode='Existing';
otherwise;
end; % switch runmode;
clear kmen;


%--- PARTIAL OR TOTAL RING WIDTH 

kmen = menu('Choose Ring-Width Type',...
    'Total Ring Width',...
    'Earlywood/Latewood Width');
switch kmen;
case 1; 
    widthmode='Total';
case 2;
    widthmode='EW/LW';
otherwise;
end;
clear kmen;


%--- PROGRAM SETTINGS (CAN HARD CODE THESE AS NEEDED

prompt={'Enter your initials','Enter the units for measurments as stored in files',...
        'Enter maximum expected length of any series (yr)'};
def={'dmm','mm x 100','1000'};
dlgTitle='Program Settings';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
person=answer{1};
units=answer{2};
maxlen=str2num(answer{3});
clear prompt def dlgTitle lineNo answer;


%--- MEASUREMENT MODE

switch runmode;
case 'Edit';
    error('Edit mode not yet implemented');
case 'Measure';
    
    kmen = menu('Choose','Start new site',...
        'Operate on existing site');
    switch kmen;
    case 1; 
        filemode='New';
    case 2;
        filemode='Existing';
    otherwise;
    end;
    
    
       
    %--- IF NEW FILE ALLOCATE
    
    switch filemode;
    case 'Existing';
        
        %-- Prompt for  name of existing rw storage file
        [file1,path1]=uigetfile('*.mat','Input file storing ring-width data');
        pf1=[path1,file1];
        if ~exist(pf1,'file'); % check that file exists
            error([pf1 ' does not exist']);
        else; % File exists; check that it has the required data
            eval(['s=load(''' pf1 ''');']);
            if ~all([isfield(s,'vlist')  isfield(s,'XEL')  isfield(s,'XT')]);
                error([pf1 ' does not contain XEL, XT and vlist']);
            end
            clear s;
            % Load the ring-width data
            eval(['load ' pf1 '  XT XEL vlist;']); 
        end; % if ~exist(pf1,'file');
        
        %-- This is an existing file'  Can add a new series or measure/edit an existing series
        kmen = menu('Choose',...
            'Start measuring a new series',...
            'Continue measuring an existing series');
        if kmen==1; % start measuring a new series in an existing storage file
            kseries = 'New';
            % Compute series index jser (one above highest existing)
            switch widthmode;
            case 'Total';
                ncurr=length(XT.id);
            case 'EW/LW';
                ncurr=length(XEL.id);
            otherwise;
            end;
            jser=ncurr+1;
            clear ncurr;
            
            % Call subfcn01.m to  name series
            nameser = subfcn01;
            
        else;  % resume measuring an existing series -- or remeasure parts
            kseries='Existing';
            % Get series id from menu
            switch widthmode;
            case 'Total';
                kmen5=menu('Choose series',XT.id);
                jser=kmen5;
                nameser=XT.id{jser};
            case 'EW/LW';
                kmen5=menu('Choose series',XEL.id);
                jser=kmen5;
                nameser=XT.id{jser};
            otherwise;
            end;
            clear kmen5;
            
        end; %  if kmen==1; % start measureing a new series in an existing storage file
        
        % Prompt for first year of measurement
        prompt={'Enter year of first value to be measured:'};
        switch kseries;
        case 'New'; % new series
            def={'1'};
        case 'Existing'; % existing series
            switch widthmode;
            case 'Total';
                nameser = XT.id{jser};
                xthis=XT.data{jser};
                yrthis=xthis(:,1);
                xthis=xthis(:,2);
                yrorigin=yrin(1);
                yrsp = yrorigin+maxlen-1;
                igo = min(yrthis)-yrorigin+1;
                isp = max(yrthis)-yrorigin+1;
                nsize = yrsp-yrorigin+1; 
                xin = repmat(NaN,nsize,1);
                xin(igo:isp)=xthis;
                yrin=(yrorigin:yrsp)';
                yrstart=max(yrthis)+1;
            case 'EW/LW';
                nameser = XEL.id{jser};
                xthis=XEL.data{jser};
                yrthis=xthis(:,1);
                xthis=xthis(:,2:4);
                yrorigin=yrin(1);
                yrsp = yrorigin+maxlen-1;
                igo = min(yrthis)-yrorigin+1;
                isp = max(yrthis)-yrorigin+1;
                nsize = yrsp-yrorigin+1; 
                xin = repmat(NaN,nsize,3);
                xin(igo:isp,:)=xthis;
                yrin=(yrorigin:yrsp)';
                yrstart=max(yrthis)+1;
            otherwise;
            end;
            
            
        otherwise;
        end;
        dlgTitle='Input Year';
        lineNo=1;
        answer=inputdlg(prompt,dlgTitle,lineNo,def);
        yrstart = str2num(answer{1});
        yrorigin=yrstart;
        clear prompt def dlgTitle lineNo answer;
                
        yrsp = yrorigin+maxlen-1;
        yrin=(yrorigin:yrsp)';
        nsize = yrsp-yrorigin+1; % allocated size of vector or matrix to store curremt series
        
        clear kmen; 
        
    case 'New'; % start new file
        
        %-- Prompt for  name of new storage file
        [file1,path1]=uiputfile('*.mat','New input file to store ring widths in ');
        pf1=[path1,file1];
        if exist(pf1,'file'); % check that file exists
            error([pf1 ' already exists']);
        else; % File exists; check that it has the required data
        end; % if ~exist(pf1,'file');
        
        % Initialize variables
        XT.id{1}=[]; % sample id -- same for total, EWW or LWW, but empty if does not exist 
        XEL.id{1}=[];
        XT.data{1}=[];
        XEL.data{1}=[];
        XT.span{1}=[];
        XEL.span{1}=[];
        XT.who{1}=[];
        XEL.who{1}=[];
        XT.when{1}=[];
        XEL.when{1}=[];
        XEL.units=units;
        XEL.summary=[];
        vlist = 'Contents of .mat file';
        vlist =char(vlist,'');
        vlist=char(vlist,'X, XT,XEL - structure variables for general inf, total ring width and earlywood/latewood width:');
        vlist =char(vlist,'   X.units -- units of measurement');
        vlist =char(vlist,'   X.summary -- series id,  first and last year of total, EWW, LWW for each each series (string matrix)');
        vlist =char(vlist,'   NEXT VARIABLES ARE CELLS, ONE ELEMENT PER SERIES, for XT, XEL');
        vlist =char(vlist,'     {id -- name of the series (e.g., PDF02A)}');
        vlist =char(vlist,'     {data -- the measurements}');
        vlist =char(vlist,'     {span -- first and last year of measurements}');
        vlist =char(vlist,'     {who -- the initials of measurer}');
        vlist =char(vlist,'     {when -- date measurements last modified}');
        
        % This is a new file; call subfcn01.m to input name of series
        jser=1;
        nameser=subfcn01;
                
        % Prompt for first year of measurement
        prompt={'Enter year of first value to be measured:'};
        def={'1'};
        dlgTitle='Input Year';
        lineNo=1;
        answer=inputdlg(prompt,dlgTitle,lineNo,def);
        yrstart = str2num(answer{1});
        yrorigin=yrstart;
        clear prompt def dlgTitle lineNo answer;
                
        yrsp = yrorigin+maxlen-1;
        yrin=(yrorigin:yrsp)';
        nsize = yrsp-yrorigin+1; % allocated size of vector or matrix to store curremt series
        switch widthmode;
        case 'Total';
            xin=repmat(NaN,nsize,1);
            % Measure
            %??[xout,yrout]=meastot(x,yrorigin,yrstart); % xout is > x 2 matrix
            xout=rand(50,1);
            yrout=(1901:1950)';
                        
            % Store measurements, span, etc
            ddate=date; %  compute current date
            XT.id{jser} =nameser;
            XT.data{jser}=[yrout xout];
            XT.who{jser} = person;
            XT.when{jser}=ddate;
            XT.span{jser}=[min(yrout) max(yrout)];
          
        case 'EL/LW';
            xin=repmat(NaN,nsize,3);
            % Measure
            %[xout,yrout]=measpart(xin,yrorigin,yrstart); % xout is ? x 3 matrix
            xout=rand(50,3);
            yrout=(1901:1950)';
                        
            ddate=date; %  compute current date
            
            % Store measurements, span, etc
                       
            XEL.id{jser} =nameser;
            XEL.data{jser}=xout(:,[1 3]);
            XEL.who{jser} = person;
            XEL.when{jser}=ddate;
            
            % Compute spans for Ew and Lw; store spans
            iE = find(~isnan(xout(:,1)));
            iL = find(~isnan(xout(:,2)));
            iC = find(~isnan(xout(:,3)));
            XEL.span{jser}=[yrout(min(iE))  yrout(max(iE))  yrout(min(iL))  yrout(max(iL))   yrout(min(iC))  yrout(max(iC)) ];
            
            % Feedback
            strmsg={['You have worked on series ' nameser],...
                    ['The measured series now covers ' int2str(XEL.span{jser})]};
            uiwait(msgbox('strmsg','Message','modal'));
            
            % Sort the revised data so that series alphabetically ordered
            XEL = subfcn02(XEL);
            
            % Update the string-array summary
            XEL = subfcn03(XEL);
            
        otherwise; % neither total nor EW/LW
        end; % switch widthmode;
        
        % Save results
        
        switch filemode;
        case 'Existing';
            kmen=menu('Choose',...
                ['Save revised data in ' pf1],...
                'Save revised data in a new file',...
                ['Abort -- original ' pf1 ' will be unchanged']);
            if kmen==1;
                eval(['save ' pf1 ' set1 -append;']);
            elseif kmen==2;
                [file2,path2]=uiputfile('*.mat','Output storage file');
                pf2=[path2 file2];
                eval(['save ' pf2  set1 ' -append;']);
            elseif kmen==3;
                disp('ABORTED ');
                return;
            end;
        case 'New';
             kmen=menu('Choose',...
                ['Save the data in new storage file ' pf1],...
                ['Abort ']);
            if kmen==1;
                eval(['save ' pf1  set1 ';']);
            elseif kmen==2;
                disp('ABORTED ');
                return;
            end;
            
            
        otherwise; % switch filemode -- saving file
        end;% switch filemode -- saving file
    otherwise;
    end;  % switch filemode (New vs Existing));
    
    %-- Sort data alphabetically -- subfcn02.m
    
    %-- Update the summary
    
    % -- Update the storage file
    
    % -- Optionally create .rw or .rwl files of revised data
    
end; %  switch runmode;


% ******  SUBFUCTION 1  -- getting and checking the series name
 
function nameser = subfcn01;

prompt={'Enter the series ID:'};
def={'XXX01A'};
dlgTitle='Series ID (site-code + tree# + core=letter (e.g., PDF01A)';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
d1=answer{1};
d2=deblank(d1);
d3=fliplr(deblank(fliplr(d2)));
nameser=upper(d3);
clear d1 d2 d3 def dlgTitle lineNo answer;

% Length of name
nlen = length(nameser);
if nlen>8;
    error(['Series ID ' nameser ' has more than 8 chars']);
elseif nlen<5;;
    error(['Series ID ' nameser ' needs at least 5 chars for a site code, tree no and core letter']);
    
end;


% First 3 chars must be letters
d1 = nameser(1:3);
if ~all(isletter(d1));
    error([nameser ' first 3 chars must be letters']);
end;
sitecode=d1;
clear d1;

% Check that have a core letter, following a tree number
d1=nameser;
d1(1:3) = []; % strip off site code
i1 = find(isletter(d1)); % find remaining letters
if isempty(i1);
    error([nameser ' needs a tree number following site code']);
else;
    treeno=d1(1:i1(1)-1);
    if length(treeno)==1;
        treeno=['0' treeno];
    elseif length(treeno)>3;
        error([nameser ' tree number should be max of 999 and max of 3 digits']);
    end;
end;

% Check core letter and for trailing core-segment number
d1 = d1(i1(1):length(d1));
i1=find(isletter(d1));
if length(d1)==1;
    coreletter=d1;
    nameser=[sitecode treeno coreletter];
else;
    % Strip the letter off, and should have only a 1-digit number
    coreletter=d1(i1(1));
    d1(1)=[];
    if length(d1)>1;
        error([nameser ' has a core-segment number after the core letter, but this should be length 1']);
    else;
        if isletter(d1); 
            error(['Last char of ' nameser ' should be a number for this case']);
        end;
        segno = d1;
        nameser=[sitecode treeno coreletter segno];
    end;
end;


% SUBFUNCTION 2 -- sort the cell data
        
function [X]=subfcn02(X)
% Sort cell storage alphabetically

% id, data , when, who, span

[X.id,j]=sort(X.id);
X.data=X.data(j);
X.when=X.when(j);
X.who=X.who(j);
X.span=X.span(j);


%--- SUBFUNCTION 3 --  update the summary
function X=subfcn03(X,k)
% update the summary

n= length(X.id); % number of series measured

if k==1; % total rw
    str1='Total Ring Width';
    str2=' ID         Total RW';
elseif k==2; % EWW/LWW
    str1='EWW/LWW';
    str2=['ID            EWW           LWW           EWW+LWW'];
end;

s1=[int2str(n) ' series  have been measured for ' str1];
s1=char(s1,'Spans of measurements are as follows');
s1=char(s1,str2);


for j = 1:n;
    spann=X.span{j};
    if k==1;;
        s2=[X.id{j} '      '  int2str(spann(1,:))];
    elseif k==2;
        s2 = [X.id{j} '     '   int2str(spann(1,:)) '     ' int2str(spann(2,:)) '     '  int2str(spann(3,:))];
    end;
    s1=char(s1,s2);
end;

X.summary=s1;

