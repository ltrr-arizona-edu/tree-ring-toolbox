function Lgood=nonan01(L)
% nonan1:  find longest stretch of ones in logical cv; 
% Lgood=nonan01(L);
% Last revised: 7-13-99
%
% A low-level function. Typical use by a calling function is to find the longest 
% unbroken string of years of climate or tree-ring variable without any missing data 
% (NaNs). Called by, among other functions, respfun0.m
%
%*** IN
%
% L (mL x 1)L    logical vector
% 
%*** OUT
%
% Lgood (mL x 1)L  longest continuous "1" stretch in L
%
%*** REFERENCES -- none
%*** UW FUNCTIONS CALLED -- none
%*** TOOLBOXES NEEDED -- none
%
%*** NOTES
%
% 

[mL,nL]=size(L);
if nL~=1;
   error('L mus be cv');
end
if ~islogical(L);
   error('L must be logical');
end


%------ Test for all OK or none aOK
if all(L);
   Lgood=L;
   return;
elseif ~any(L);
   Lgood = logical(zeros(mL,1));
else
   L1 =[0; L; 0];
   i1 = find(diff(L1)==1);
   i2 = find(diff(L1)==-1);
   i3 = i2-i1;
   [i3max,j]=max(i3);
   igo = i1(j);
   isp = igo+i3(j)-1;
   L2 = zeros(mL+1,1);
   L2(igo:isp)=1;
   L2(mL+1)=[];
   Lgood = logical(L2);
end








