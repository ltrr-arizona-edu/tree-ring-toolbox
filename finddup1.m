function [J,D]=finddup1(S,n,kopt)
% funddup1: find duplicate names in a string matrix
% [J,D]=finddup1(S,n,kopt);
% Last revised 3-4-00
%
% Find and report rows of a string matrix that are duplicates.  Control over 
% matching only first few chars, or the chars up to a '.'.
%
%*** INPUT 
% 
% S (nS x mS)s  string matrix to be checked (e.g., a list of file names)
% n (1 x 1)i  the number of characters to compare (see notes)
% kopt(1 x 1)i options
%	kopt(1)==1  match based on leftmost n chars, ignoring leading blanks
%	       ==2  match chars to left of '.'.  Ignore the "n" spec
%
%*** OUTPUT
%
% J (? x 2)i  the rows of S that are dupes; if no dupes, J is []
% D (? x ?)s  strings for rows in the dupes; if no dupes, D is empty string
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED -- NONE
%*** TOOLBOXES NEEDED -- NONE
%
%*** NOTES
%
% The strings from rows of S are left-deblanked before comparison, so that, say,
%   fre001.crn   and
%  fre001.crn   will show as matches
%
% The '.' option (kopt(1)==2) is useful if the rows of S are filenames, as you 
% compare the part before the suffix.


if ~ischar(S);
   error('S must be string matrix');
end;


% Initialize
J=[];

% Left justify S
S1=fliplr(S);
S2=cellstr(S1);
S2=deblank(S2);
S2=char(S2);
[mS2,nS2]=size(S2);
D=blanks(nS2);
if nS2<n & kopt(1)==1;
   error('Col size of S2 must be as large as n if using the n-chars option');
end;

S2=fliplr(S2);
SC=cellstr(S2);
L=zeros(mS2,1);

% Loop over rows  of S
for m = 1:mS2;
   %disp(m);
   T2=S2;
   T2(m,:)=[];
   s2=S2(m,:);
   if kopt(1)==1; % n option
      s2=s2(1:n);
   else;
      if isempty(findstr(s2,'.'));
         error([s2 ' does not have a period in it']);
         s2=strtok(s2,'.');
      end;
   end;
   j1=strmatch(s2,T2);
        
   if ~isempty(j1);
      D=char(D, s2);
      J=[J; m];
            
   end;
end
D(1,:)=[];
   