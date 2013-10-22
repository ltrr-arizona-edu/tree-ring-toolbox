function [s,t,c,f,flnm2]=treefn(flnm1)
%
% From a .RW, .EWW, or .LWW file, get the 3-letter site code,
% tree number, core letter, and file-type; and build a 
% concise file name.
%
%****************** IN ARGS *******************
%
% flnm1 (1,?)s   input file name
%
%
%******************** OUT ARGS *******************
%
% s (1 x 3)s   3-char site code
% t (1 x 1)i   tree number
% c (1 x ?)s   core letter (possibly with a number on end)
% f (1 x 3)s   data type -- total ringwidth, earlywood width or
%				latewood width
% flnm2(1 x ?)s  output file name
%
%******************* NOTES ******************************
%
% Why?  File naming convention for .RW, .EWW, .LWW files is so
% loose that couldn't easily test whether I got the right .RW file
% to match, say, the EWW and LWW files for the same core.  
% treefn.m extracts the meat from the file names, and builds a
% meaty filename.  For example, trees 020, 02 , 2 all become
% tree 2.
%
% Three-letter site code is assumed to be first 3 chars of filename
% Next numeric part is tree number
% Next part, up to the "." is the core id
% Suffix indicates the file type 
%  EWW= earlywood width
%	LWW= latewood width
%	RW or TRW == total ring width
%
% "xe" or "xl" part of core id is disregarded 
%
%
% Examples:
%
%     flnm1            s    t  c    f        flnm2
%
%    pdf101a.RW      pdf  101 a    TRW --> pdf101a.TRW
%    pdf01axe.EWW    pdf  1   a    EWW --> pdf1a.EWW
%	  pdf1a.RW        pdf  1   a    TRW --> pdf1a.TRW
%

x=flnm1;
n1=length(x);



%*********** site code

s = x(1:3);


%**************  suffix

sfx = x((find(x=='.')+1):n1);
if length(sfx)==2
	if   all(sfx=='RW')
		f = 'TRW';
	else
		error('Only allowable 2-char suffix is RW')
	end
elseif length(sfx)==3;
	if  all(sfx=='TRW')
		f='TRW';
	elseif all(sfx=='EWW')
		f='EWW';
	elseif all(sfx=='LWW')
		f='LWW';
	else
		error('Only allowable 3-char suffixes: TRW,EWW,LWW')
	end
else
	error('Suffix must be of length 2 or 3')
end

%************** core id

x = strtok(x,'.');  % part before the '.'
x(1:3)=[];  % cut off first 3 chars

n3=length(x);

f1 = find(isletter(x));  % first letter should be start of core indicator
f1=f1(1);

c = x(f1:length(x));
nc = length(c);
if nc>1
	f3=find(c=='X');
	if ~isempty(f3)
		if(f3(1)>1)
			c(f3(1):nc)=[]; % lop off what will be 'XE' or 'XL'
		end
	end
end


t =str2num(x(1:(f1-1)));
t = num2str(t);

flnm2=[s t c '.' f];
