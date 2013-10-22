% mdncull.m   Given an array with one column holding a type indicator and
%  another column a variable, computes the minium, median, and maximum by
%  type.

%******************  PRELOADS  *******************************************

% 1.  The data array, by any name.  Will be screen prompted to rename
%     as X within pgm


%*******  Data Setup

clear X
f=input('Filename, without extension: ','s');
X=eval(f);  % data array in X

p1=input('INDICATOR COLUMN: ');
p2=input('VARIABLE COLUMN: ');

kk=1;  % toggle for while loop

while kk==1
	k=menu('CHOOSE A SPECIES','YES','QUIT');
	if k==1  % choose another species
		i=input('SPECIES #: ')
		clear  x y;
		
		I=X(:,p1)~=i;
		y=X(:,p2);     % the variable
		y(I)=[];     % delete all rows not corresp to desired vbl
		clc;
		disp('INDICATOR  SAMPLSIZE   MINIMUM   MEDIAN    MAXIMUM')
		disp('');
		disp([i  length(y)  min(y)  median(y)  max(y)]);
	else
		kk=0;   % end loop, no more work
	end

end


