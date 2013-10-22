function [zname,nsers]=namenum(nmline)
%
% string matrix of core or tree names with sequential numbers
%
%
% D Meko 8-10-95
%
%
%************************** IN ARGS *************************
%
% nmline (1 x ?)s  line of names (8-or-fewer characters)
%		separated by white space

%
%************************* OUT ARGS ************************
%
% zname (? x 11)s matrix of names with sequence numbers at
%		at right
% nsers (1 x 1)i number of series in matrix
%
%
%



nmline=[' ' nmline '        '];
l=length(nmline);

% Find col position of first char in each name
L1 = ~isspace(nmline);
d1 = diff(L1);
d1  = [' ' d1];
L2 = d1 ==1;
I2 = find(L2);  % index of first char

% Find col position of last char in each name
L3 = d1==-1;
I3=find(L3);

% If first "name" begins with a number, it is not really a
% name but a "number of columns" indicator as put on the line
% by program "YUX" in the program library.  If such an
% indicator exists, replace it with blanks and adjust I2 and I3
% accordingly
if ~isletter(nmline(I2(1)))
	nsers=sum(L2)-1;
	I2(1)=[];
	I3(1)=[];
else
	nsers = sum(L2);  % number of "first characters"
end


vblank=(blanks(nsers))';
zname1 = vblank(:,ones(8,1));
for i = 1:nsers	
	nchars = I3(i)-I2(i)+1;
	zname1(i,1:nchars) = nmline(I2(i):I3(i));
end




%*************  Build the Sequence No. Part of the Matrix
%
s1a = '1234567890';
s1b = s1a(ones(nsers,1),:);
s1c = (s1b)';
[m1,n1]=size(s1c);
S1 = reshape(s1c,(m1*n1),1);
S1 = S1(1:nsers);

s2a = '0123456789';
s2b = s2a(ones(10,1), :);
[m2,n2]=size(s2b);
s2c = reshape(s2b,(m2*n2),1);
S2 = s2c(2:nsers+1);

S3 = (blanks(nsers))';


S = [S3 S2 S1];


%*************************** Splice names and seq numbers

zname = [zname1 S];
