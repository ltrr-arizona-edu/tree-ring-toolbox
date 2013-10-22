function T=strip1(S,txt);
% strip1: strip specified text string off end of each row of matrix S and replace with blanks
% CALL: T=strip1(S,txt)
%
%************** IN *******************
%
% S (mS x nS)ch   character matrix of, say, core names
% txt (1 x ?)ch  character string to replace with blanks
%
%***************** OUT ****************
%
% T (mT x nT)ch   Like S, but with ending chars equal to 'txt' replaced 
%
%******************* NOTE **************
%
% Misc function written to remove 'XE' from tail end of earlywood core names so
% that functions in standardization and tree identification worked correctly

[mS,nS]=size(S);
if ~ischar(S);
   error('S must be char ');
end

% Initialize T as S
T=S;

for n = 1:mS;
   c=strtok(S(n,:));
   if length(c)>=length(txt);
      isp = length(c);
      igo = length(c)-length(txt)+1;
      if strcmp(txt,c(igo:isp)); % if the ending chars match txt
         T(n,igo:isp)=blanks(length(txt));
      else
      end
   end
end
