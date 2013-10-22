function L=menudm1(tit,names,nmax)
% menudm1: interactive selection of subset of menu items
% CALL:  I=menudm1(names,nmax);
%
% Meko 3-17-98
%
%*****************  IN
%
% names {}s    1x? cell matrix of names to be picked from
%
%****************** OUT
%
% L (1 x ?)L  pointer to selected names (1) or not selected (0)
% nmax (1 x 1)i  maximum number of selected names allowed


ntot = size(names,2); % number of items to pick from
if nmax>ntot;
   error('nmax bigger than ntot');
end

L=logical(zeros(ntot,1));

for n = 1:ntot;
   names{n}=['N-' names{n}];
end
names{ntot+1}='Satisfied';

npick = 0;
kwhile1=1;
while kwhile1==1;
   kmen1 = menu(tit,names);
   if kmen1~=ntot+1;
      if L(kmen1);
         names{kmen1}(1:2)='N-';
         L(kmen1)=0;
      else
         names{kmen1}(1:2)='Y-';
         L(kmen1)=1;
      end
   else
      if sum(L)>nmax;
         msgbox('You picked one too many series. Unmark one!');
      else
         kwhile1=0;
      end
      
   end
   if sum(L)>nmax;
      msgbox('You picked one too many series. Unmark one!');
   end
   
end

      