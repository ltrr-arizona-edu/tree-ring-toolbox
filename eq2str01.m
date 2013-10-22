function str = eq2str01(th)
% eq2str01:  ARMA model equation to string; for displaying equation
% CALL str=eq2str01(th);
%
% Meko 2-20-97
%
%************ IN 
%
% th -- "theta" format model info returned by functions such as ar.m
%
%************ OUT 
%
% str (1 x ?)s  string with equation



%************  FIND OUT AR ORDER AND MA ORDER

na = th(1,4); % AR order
nc = th(1,5); % MA order

% Make sure nu zero
ntemp = th(1,[3 6 7]);
if ~all(ntemp==0);
   error('Row 1 of th matrix implies not an ARMA, AR, or MA model');
end


%***************  BUILD left-hand string
str1='\ity_{t} '; % initialize left hand side
if na>0;
   ca = th(3,1:na);
   str1a=' ';
   for n = 1:na;
      c =ca(n); % nth ar coef
      str1a1 = sprintf('%4.3f',c);
      kterm = int2str(n);
      str1a2 = ['\ity_{\itt-\rm' kterm '}'];
      if c>=0;
         sgn = '+';
      else
         sgn = ' ';
      end
      str1a=[str1a sgn str1a1 str1a2];
   end
   str1=[str1 str1a ' = '];
else
   str1=[str1 '= '];
end



%*************** BUILD right-hand string
str2='\ite_{t} ';  % Will always have the e(t) term

if nc>0;
   cc = th(3,(na+1):(na+nc));
   str2a = ' ';
   for n = 1:nc;
      c = cc(n); % nth ma coef
      str2a1 = sprintf('%4.3f',c);
      kterm = int2str(n);
      str2a2 = ['\ite_{\itt-\rm' kterm '}'];
      if c>0;
         sgn='+';
      else
         sgn=' ';
      end
      str2a=[str2a sgn str2a1 str2a2];
   end
   str2 = [str2 str2a];
else
   str2 = [str2];
end


str=[str1 str2];



