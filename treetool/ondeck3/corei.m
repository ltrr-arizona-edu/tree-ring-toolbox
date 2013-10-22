function corei
% corei:  compute standard and AR-residual core tree-ring index
% corei;
% Last revised: 9-2-99
%
% Computes core index from core ring width and previously fit 
% growth trend.  Also AR-models the ring-width series and computes 
% residual core index. Second in sequence of functions, following crvfit.m, which 
% is used to fit the growth trend.
%
%*** INPUT
%
% User prompted to point to the .mat storage file which has the ring-width 
% and curve-fit information.
%
%*** OUTPUT
%
% User prompted to point to .mat storage file to store the computed indices in.  
% Typically, this is the same as the input .mat file with the ring width. In 
% this case, the computed indices are appended to the .mat file.
%
%*** REFERENCES -- NONE
%
%*** UW FUNCTIONS CALLED
% 
% armod
%   resid1
%   acf
%   pacf
% stats1
%   meansen1
%   basic1
%   rtree1
%     treenum
%     maskind
%   rtree2




%*** TOOLBOXES NEEDED
%
%*** NOTES

clear all

% Prompt for the .mat file name

% Get the .mat filename interactively; load file
[flmat,path1]=uigetfile('*.mat','Input .MAT filename ?');
pf1=[path1 flmat];
flold=flmat; % Will need this file name later
eval(['load ',pf1]);
clear flmat;


[ns,dum1]=size(nms);  % Store the size of nms string matrix

if ~exist('S'),
  error('No curves have been test-fitted yet');
end
if all(S(:,1)==zeros(ns,1)),
  error('First column of S is all zero');
end

IX=ones(length(X),1)*NaN; % Standard indices will be stored in same-size 
  % vector as ring widths

nblk=5; % Maximum possible number of blocked intervals
scid=1;
disp('Growth-curve fitting');

while scid<=ns,
     if S(scid,1)~=0,
          
     % Get the Full-length ring-width series and its years vector 
     yrv=(yrs(scid,1):yrs(scid,2))';
     xv=X(yrs(scid,3):yrs(scid,3)+yrs(scid,2)-yrs(scid,1));
     
     % Get the fitted growth curve
     gv=G(yrs(scid,3):yrs(scid,3)+yrs(scid,2)-yrs(scid,1));
     
     % Get row index within yrv (and xv and gv) of beginning and ending years
     % of fit-segment as marked in cfit1.m and stored in S
     xind1=find(yrv==S(scid,10));
     xind2=find(yrv==S(scid,11));
     
     % Get the coresponding years, ring-width data, and growth trend data
     yrvn=S(scid,10):S(scid,11);
     xvn=xv(xind1:xind2);
     gvn=gv(xind1:xind2);
     
     % Make a logical pointer to the non-NaN elements in gvn and xvn and check 
     % that the same elements are non-NaN in both
     L1=~isnan(xvn);
     L2=~isnan(gvn);
     if any(L1 & ~L2) | any(L2 & ~L1); 
        error('All non-NaN elements of gvn and xvn are not in same rows');
     else;
        L=L1;
        clear L1 L2;
     end;
     
     % Check for zero values in growth curve -- problems for division
     L1 = gvn==0;
     if any(L1);
        error(['growth curve for series ' int2str(scid) ' has zero values']);
     end;
     
     %--- Comppute the index by the ratio method
     
     rwin1=xvn./ gvn;
     
     % Adjust ring-width index by shifting to mean 1.0
     mnx = nanmean(rwin1);
     rwin2 = rwin1 -mnx + 1.0;
     
     % Store index
     IX(S(scid,12):S(scid,12)+S(scid,11)-S(scid,10))=rwin2;
     disp(['   '   nms(scid,:) ' fit']);

  end
  scid=scid+1;
end		% End of scid while loop


%******************** AR MODELING AND GENERATION OF RESIDUAL CORE INDICES

disp(' ');
disp('AR Modeling');

prompt={'Maximum Order:'};
def={'2'};
dlgTitle='Maximum Allowable AR Model Order';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
nnmx=str2num(answer{1});

disp(['Maximum AR order to consider in prewhitening model: ' int2str(nnmx)]);

[EX,P]=armod(IX,S,yrs,nms,nnmx);

%******************* STATISTICS OF CORE INDICES *********************

[CE1,CE2,CE3,CS1,CS2,CS3]=stats1(IX,EX,yrs,nms);


%*********************************


% Save the vectors in a .mat file
setnew=' IX EX P CE1 CE2 CE3 CS1 CS2 CS3 '; % New variables to be saved
nsv=menu('Save variables?','Add New Variables to original file?',...
   'Store new variables only in a new file?','QUIT');
if nsv==1,
  eval(['save ' pf1 setnew ' -append']);
elseif nsv==2,
   [ofmat,path2]=uiputfile('*.mat','New .MAT filename to store new stuff in: ');
   pf2=[path2 ofmat];
  eval(['save ',pf2, setnew]);
end


% End of file
