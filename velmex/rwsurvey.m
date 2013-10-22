function rwsurvey(X,Y)
% rwsurvey:  build .txt files summaring time coverage of series
%  rwsurvey(X,Y);
% Last modified 5-28-03
%
% Build .txt files summaring time coverage of series; function called by rwmeas
%
%*** INPUT 
%
% X: structure with measured-total-width data (see rwmeas)
% Y: structure with partial width data
% rwlset.  structure of .rwl file sets, with fields:
%
%*** OUTPUT
% 
% No arguments
% Optional output .txt files with series names and time coverage
%
% 
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% List of series can be sorted 
% 1) alphabetically by id
% 2) by start year
% 3) by length of series


txtmen1 = {'Measured Total ','Partial','Quit'};
txtmen3={'Alphabetical by id','Numerical by start year','Numerical by end year','Numerical by series length'};
txtmen2={'Earlywood','Latewood','Computed Total'};

kwh1=1;
while kwh1;
    kmen1=menu('Choose type of data to survey',txtmen1);
    if kmen1==3;
        kwh1=0;
    elseif kmen1==1; % survey measured total
        W=X;
        datatype='Total';
        k1=1;
        k2=[];
        kmen3=menu('Choose sort criterion',txtmen3);
        txthead = [txtmen1{kmen1}  '/' txtmen3{kmen3}];
        [kflag1,id,D]=subfcn01(W,k1,k2);
    elseif kmen1==2; % survey partial or computed total
        W=Y;
        datatype='Partial';
        k1=2;
        kmen2=menu('Choose data',txtmen2);
        kmen3=menu('Choose sort criterion',txtmen3);
        txthead = [txtmen1{kmen1} '/' txtmen2{kmen2} '/' txtmen3{kmen3}];
        k2=kmen2;
        [kflag1,id,D]=subfcn01(W,k1,k2);
    end;
    
    if kflag1==0;
       
        % no series
    else;
        
        if kmen1~=3; 
            
            % Get sort row index
            if kmen3==1; % alpahabet
                [S,js]=sort(id);
            elseif kmen3==2; % start yr
                [S,js]=sort(D(:,1));
            elseif kmen3==3; % end yr
                [S,js]=sort(D(:,2));
            else; % length
                [S,js]=sort(D(:,3));
            end;
            
            % Sort;
            A=char(id(js));
            B=D(js,:);
            [mA,nA]=size(A);
            bb=bnkcols(mA);
            C=[A bb.n2  num2str(B)];
            
            head1='ID       First Last    N';
            C1=char(txthead);
            C1=char(C1,['Made by rwmeas/rwsurvey on ' datestr(date)]);
            C1=char(C1,blanks(3));
            C1=char(C1,head1);
            C1=char(C1,blanks(3));
            C=char(C1,C);
            
            
            % Write tables to file
            
            strmess1={'You will be prompted for a name of an output file.',...
                    'You can open that file with vedit.',...
                    'Do not give the file a name that is the same as a ',...
                    '.txt file currently open in vedit.  Good strategy is to ',...
                    'use the same name (e.g., aaasur.txt) over and over',...
                    'overwriting each time',' ',...
                    'If you get error message:   ??? Error using ==> fprintf , Invalid fid.',...
                    'you forgot to close an old version of the same file in vedit'};
              
             uiwait(msgbox(strmess1,'Message','modal'));
            
            [file1,path1]=uiputfile('aaa*.txt','Name for output .txt file');
            pf1=[path1 file1];
            
            fid1=fopen(pf1,'w');
            
            [mC,nC]=size(C);
            eval(['fmt1=''%' int2str(nC) 's\n'';'])
            for n =1:mC;
                s = C(n,:);
                fprintf(fid1,fmt1,s);
            end;
            fclose(fid1);
            uiwait(msgbox({[pf1 ' successfully written!'],...
                    'You may now request additional files'},'Message','modal'));
            
        end;
        
    end;
end; % while kwh1;



function [kflag1,id,D]=subfcn01(S,k1,k2)
if k1==1; % measured total
    nid = length(S.id);
    if nid==1 & isempty(S.id{1});
        D=[];
        id=[];
        kflag1=0; % no series
        uiwait(msgbox(['No measured-total series'],'Message','modal'));
        return;
    else; % at least one series
        kflag1=1;
        id = S.id;
        A=repmat(NaN,nid,3); % to hold first year, last year, number of years
        for m = 1:nid;
            a = S.span{m};
            b= a(2)-a(1)+1;
            A(m,:)=[a b];
        end;
    end;
elseif k1==2; % partial 
    nid = length(S.id);
    if nid==1 & isempty(S.id{1});
        D=[];
        id=[];
        kflag1=0; % no series
        uiwait(msgbox(['No partial series'],'Message','modal'));
        return;
    else; % at least one series
        kflag1=1;
        id = S.id;
        A=repmat(NaN,nid,3); % to hold first year, last year, number of years
        for m = 1:nid;
            a = S.span{m}(k2,:);
            b= a(2)-a(1)+1;
            A(m,:)=[a b];
        end;
    end;
end;
D=A;


