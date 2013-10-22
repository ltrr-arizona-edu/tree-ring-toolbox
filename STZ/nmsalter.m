function nmsalter
% nmsalter:  Alters nms from a .mat tree-ring storage file to get rid of strange forms
% nmsalter;
% Last revised 9-28-99
%
% In running the crvfit series of standardization programs, it is assumed that only 
% one letter will be after the core number (e.g., DCY12A).  Forms such as 
% DCY12AXL  are not acceptable, and will cause problems later in functions that 
% attempt to decipher tree number and core number from the code.  Nmsalter.m 
% circumvents this problem.  Run this (once only required) after you convert the 
% .rwl file to a .mat file with rwlinp.m.
%
%*** IN
%
% You point to a .mat storage file
%
%*** OUT
%
% The variable nms is revised to get rid of any letter that follow the core letter.
%
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES

%-- Get the input .mat file; nms is all that is needed
[file1,path1]=uigetfile('*.mat','Input .mat file holding nms');
pf1=[path1 file1];
eval(['load ' pf1 ' nms;']);

% logical pointer to letters
L=isletter(nms);
ncore = size(L,1);
nlength=size(L,2);

% Check that first 3 chars are letter
L1=any(any(~isletter(nms(:,1:3))));
if L1;
   error('nms must have all letters in first 3 cols of each name');
end;

% Check that 4th char is not a letter
nm4= nms(:,4);
if   any(isletter(nm4));
   error('4th character in some row of nms is a letter; should be a number');
end;

% Loop over the names.  
% Find the col that ends the tree number,...
% Check Check that in no name is there a letter following the 
% tree number.  Sufficient to check that 
for n =1:ncore;
   nm=nms(n,:);
   nmcode = nm(1:3);
   
   % Pull the segment that begins with the tree number
   nm1 = nm(4:nlength);
   iletter = find(isletter(nm1));  % pointer to positions of letters
   i1=iletter(1); % position of first letter
   nmstart= [nmcode nm1(1:i1)];  % first part of name -- thru core letter
   nstart=length(nmstart);  % length of first part of name
   
   % Assign the second part of the name -- the part after the legit core letter
   if nstart==nlength; % nothing left to do
   else;
      nmrest = nm((nstart+1):nlength); % rest of name
      nrest=length(nmrest); % how many chars in second part
      
      % Replace any letters with blanks 
      L1=isletter(nmrest);
      numlet = sum(L1);
      if numlet>0; % illegal letters, must replace
         nmrest(L1)=[];
         nmrest=[nmrest blanks(numlet)];
      else;
      end;
      nmnew = [nmstart nmrest];
      nms(n,:)=nmnew;
         
   end;
end;


% Put the revised nms back in the .mat storage file
eval(['save ' pf1 ' nms -append;']);
          
      
      
      
      
   



