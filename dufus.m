function [X,strout]=respfn0a
%respfn0a:  subfunction of respfun0 that gets the tree-ring chronology

kmen1 = menu('Type of Source File for Tree-Ring Series',...
   'ITRDB-style .crn', ...
   '2-column ascii, with year in col 1',...
   'Matlab chrononology file');
kchron = menu('Chronology type',...
   'Standard',...
   'Residual',...
   'ARSTAN');
switch kchron;
case 1;
   chrontyp='Standard Chronology';
case 2;
   chrontyp ='Residual Chronology';
case 3;
   chrontyp = 'ARSTAN Chronology';
end

%-------- Get and store tree-ring data
switch kmen1;
case 1; % .crn file
   [file1,path1]=uigetfile('*.crn','Input .crn file with tree-ring chron');
   pf1 = [path1 file1];
   [x,s,yr]=crn2vec2(pf1); % call user-written function to get the data
   % Lop-off pre-1850 data because not needed for response function, generally
   Ltemp = yr<1850;
   if any(Ltemp);
      x(Ltemp)=[]; yr(Ltemp)=[];
   end
   % Store the data
   X=[yr x];
   % Cleanup
   clear x s yr Ltemp;
case 2; % 2-column ascii file
   [file1,path1]=uigetfile('*.dat','Input .dat file with tree-ring chron');
   pf1 = [path1 file1];
   eval(['X=load(pf1);']);
   Ltemp = X(:,1)<1850;
   if any(Ltemp);
      X(Ltemp,:)=[];
   end
case 3;  % chronology .mat storage file
   [file1,path1]=uigetfile('*.mat','Input .mat file with tree-ring chron');
   pf1 = [path1 file1];
   % Chronologt will be first col in ZI or ZE, depending on if standard or
   % residual chron, and year will be in yrZI or yrZE
   if kchron==1; % standard chron
      eval(['load ' pf1 ' ZI yrZI;']);
      yr = yrZI;
      x = ZI(:,1);
      clear ZI yrZI;
   elseif kchron==2; % residual chron
      eval(['load ' pf1 ' ZE yrZE;']);
      yr = yrZE;
      x = ZE(:,1);
      clear ZE yrZE;
   else;
      error('kchron must be 1 or 2 for chrons in .mat storage files');
   end
   Ltemp = yr<1850;
   if any(Ltemp);
      x(Ltemp)=[]; yr(Ltemp)=[];
   end
   X = [yr x];
        
end

strout='TREE-RING SERIES';
strout = char(strout,['  File = ' pf1]);
strout=char(strout,['  Version = ' chrontyp]);
