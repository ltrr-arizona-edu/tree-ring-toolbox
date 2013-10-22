% script file addtxt
% addtxt:  add or replace a line in a variable-definition matrix varlist
% addtxt
% Last revised 5-4-99
%
% Utility script that facilitates updating a variable-definition list.
% Assumes the user is storing the definitions of variables in a char matrix
% called vlist.  Typical scenario is that variables are stored in a .mat file along
% with a corresponding vlist. 

if exist('kappend')==1 | exist('texttoadd')==1 | exist('mrowsize')==1 | ...
       exist('rowtoreplace')==1;
   error('One or more contradictory names already in workspace');
end

if exist('vlist')~=1;
   error('vlist does not exist in workspace');
end

mrowsize=size(vlist,1); % number of rows in varlist
vlistnew=vlist;

%-- get new text
texttoadd = input('text: ','s');


%-- append or substitute
kappend=questdlg('Append to vlist?');
switch kappend;
case 'Yes';
   vlistnew=char(vlistnew,texttoadd);
   vlist = vlistnew;
otherwise;
   rowtoreplace=input('Row to replace: ');
   if rowtoreplace>mrowsize;
      error(['vlist has fewer than ' int2str(rowtoreplace) ' rows']);
   end
   disp(vlist(rowtoreplace,:));
   disp(texttoadd);
   
   
   kgoforit=questdlg('Is that OK?');
   switch kgoforit;
   case 'Yes';
      if rowtoreplace==1;
         if mrowsize==1;
            vlistnew=texttoadd;
         else;
            vlistnew=char(texttoadd,vlistnew);
         end
         vlist=vlistnew;
      else;
         vlistnew=char(vlist(1:(rowtoreplace-1),:),texttoadd);
         vlistnew=char(vlistnew,varlist((rowtoreplace+1):mrowsize,:));
         vlist=vlistnew;
      end
   otherwise;
   end
   
   
end

clear rowtoreplace mrowsize kappend texttoadd kgoforit vlistnew
