function [Y,loggy]=rwedit(X)
% Y=rwedit(X)
%
% Edit ring-width measurements by adjusting for locally
% absent (LA) or false rings
%
% D Meko 1-1-96
%
%***************** IN ARGS 
%
% X (mX x nX)r   years (col 1) and ringwidths (col 2), before editing
%
%****************** OUT ARGS
%
% Y (mY x nY)r  years and ringwidths, after editing
% loggy (1 x 15)s string indicating what the editing was
%     Example:  'F2RF(1812)-0.25'  means 1812 had a misinterpreted
%		false ring that should have been real, and that the
%		ring has been split such that 1/4 (0.25) of the input
%		1812 has been assigned to 1812, and the remainder to 1813 in
%		the edited file.  If no editing, loggy is assigned 15 'x' 
%
%***************** UW FUNCTIONS CALLED -- NONE
%
%*************** NOTES
%
% Why? For interactive editing in rwlook.m suite.
%
% rwedit screen prompts for the type of adjustment to make. Eight
% possible adjustments can be made:
% 1 - INSL insert LA ring and adjust later dates
% 2 - INSE insert LA ring and adjust earlier dates
% 3 - DELL delete LA ring and adjust later dates
% 4 - DELE delete LA ring and adjust earlier dates
%
% 5 - convert mis-identified false ring to real and adj later dates
% 6 - convert mis-identified false ring to real and adj earlier dates
% 7 - convert mis-identified real ring to false and adj later dates
%		Sum of width for specified year t and for year t+1 in original file
%		will show up as width for year t in new file
% 8 - convert mis-identified real ring to false and adj earlier dates
%		Sum of width for specified year t and for year t+1 in original file
%		will show up as width for year t+1 in new file
%
% REVERSING OPERATIONS
% INSL(t)  -- DELL(t)
% INSE(t) --  DELE(t)
% F2RL(t) -- R2FL(t)
% F2RE(t) -- R2FE(t-1)  ***** NOTE THIS ****
%
% R2FL and R2FE cannot easily be reversed, because you would
% need to know the ratios of ringwidths in t and t+1 to recover
% using F2RL and F2RE

% Check consistency of X
% X must be 2-column matrix, with elements in first column increasing
% by increments of 1 down the column (year or nominal year), and
% elements of column 2 non-negative and less than 999.
[mX,nX]=size(X);
if nX~=2;
	error('X must be 2-col matrix')
end
d1=diff(X(:,1));
d2=d1==1;
if ~all(d2),
	error ('Column 1 of X must increment by 1')
end
L1=X(:,2)<0 | X(:,2)>=999;
if any(L1),
	error('Measurements in col 2 of X must be + and  <999')
end

% mX is number of years of ring width in X

a=NaN;

% String matrix for getting entry to edit log
elog=[...
'INSL(';'INSE(';'DELL(';'DELE(';...
'F2RL(';'F2RE(';'R2FL(';'R2FE(';'CLST('];
decf='     ' % slot for decimal fraction assigned to earliers
	% of ring pairs in R2F options

% Store year and rw measurements in vectors
yrx=X(:,1);
x=X(:,2);
yrlast=max(yrx); % last year of ring width as input

kk1 = menu('Type of Editing','Change Last Year','Other');
if kk1==2;
   % Screen prompt for key year
   kyr = input('Change value for this year: ');
   if kyr < yrx(1) | kyr> yrx(mX),
      error('Specified year outside range in X ')
   end
   if kyr==yrx(1) | kyr==yrx(mX),
      error('Specified year is first or last year of X')
   end
   
   % Compute row indices to X of parts of X not to be changed:
   % these are the parts before and after the specifed year
   I1 = yrx<kyr;
   I3 = yrx>kyr;
   x1=x(I1);
   x3=x(I3);
   xkey = x(yrx==kyr); % ringwidth for specified key year
elseif kk1==1;
end

   
   
k1=1;        % while control for level-1 menu
while k1<=9;  
   if kk1~=1;
      k1=menu('Choose one: ',...
         'INS LA & SHIFT LATER',...
         'INS LA & SHIFT EARLIER',...
         'DEL LA & SHIFT LATER',...
         'DEL LA & SHIFT EARLIER',...
         'FALSE TO REAL & SHIFT LATER',...
         'FALSE TO REAL & SHIFT EARLIER',...
         'REAL TO FALSE & SHIFT LATER',...
         'REAL TO FALSE & SHIFT EARLIER',...
         'CHANGE LAST YEAR',...
         'OOPS -- NO EDITING NEEDED');
   else
      k1=9;
   end
   
      



if k1==1;   % INSL
	yry=[yrx; yrx(mX)+1];
	y = [x1; 0; xkey; x3];
	Y=[yry y];
	loggy=[elog(1,:) int2str(kyr) ')' '     '];
	k1=10;
elseif k1==2; % INSE
	yry=[yrx(1)-1; yrx(1:mX)];
	y = [x1; xkey; 0;  x3];
	Y=[yry y];
	loggy=[elog(2,:) int2str(kyr) ')' '     '];
	k1=10;
elseif k1==3; % DELL
	if xkey~=0,
		error('Specified ringwidth not zero ')
	end
	yry=[yrx(1:mX-1)];
	y = [x1; x3];
	Y=[yry y];
	loggy=[elog(3,:) int2str(kyr) ')' '     '];
	k1=10;	
elseif k1==4; % DELE
	if xkey~=0,
		error('Specified ringwidth not zero ')
	end
	yry=[yrx(2:mX)];
	y = [x1; x3];
	Y=[yry y];
	loggy=[elog(4,:) int2str(kyr) ')' '     '];
	k1=10;	
elseif k1==5; % F2RL
	yry=[yrx; yrx(mX)+1];
	% Split the ring
	fracte=input('Decimal fraction for first year of pair: ');
	if fracte<=0 | fracte>1.0,
		error('fracte must be between 0 and 1')
	end
	xE = round(xkey*fracte); % first part
	xL = xkey-xE; % remainder
	y = [x1; xE; xL; x3];
	Y=[yry y];
	s1 = num2str(fracte);
	lens1=length(s1);
	% Pad s1 with blanks, if needed to bring to length 4
	if lens1<4,
		n1=4-lens1;
		b1=blanks(n1);
		s1=[s1 b1];
	end
	s2=[':' s1(1:4)];
	loggy=[elog(5,:) int2str(kyr) ')' s2];
	k1=10;
elseif k1==6; % F2RE
	yry=[(yrx(1)-1); yrx];
	% Split the ring
	fracte=input('Decimal fraction for first year of pair: ');
	if fracte<=0 | fracte>1.0,
		error('fracte must be between 0 and 1')
	end
	xE = round(xkey*fracte); % first part
	xL = xkey-xE; % remainder
	y = [x1; xE; xL; x3];
	Y=[yry y];
	s1 = num2str(fracte);
	lens1=length(s1);
	% Pad s1 with blanks, if needed to bring to length 4
	if lens1<4,
		n1=4-lens1;
		b1=blanks(n1);
		s1=[s1 b1];
	end
	s2=[':' s1(1:4)];
	loggy=[elog(6,:) int2str(kyr) ')' s2];
	k1=10;
elseif k1==7; % R2FL
	yry=[yrx(1:mX-1)];
	y = [x1; xkey+x3(1); x3(2:length(x3))];
	Y=[yry y];
	loggy=[elog(7,:) int2str(kyr) ')' '     '];
	k1=10;
elseif k1==8; % R2FE
	yry=[yrx(2:mX)];
	y = [x1; xkey+x3(1); x3(2:length(x3))];
	Y=[yry y];
	loggy=[elog(8,:) int2str(kyr) ')' '     '];
   k1=10
elseif k1==9; % New last year
   prompt={'Enter last year:'};
   def={num2str(yrlast)};
   tit=['Change the last year of ring width from ' int2str(yrlast)];
   lineNo=1;
   answer=inputdlg(prompt,tit,lineNo,def);
   yrnew = str2num(answer{1});
   loggy=[elog(9,:) int2str(yrnew) ')' '     '];
   Y = [((yrnew-mX+1):yrnew)' X(:,2)];
   k1=10;
elseif k1==10; % no editing needed
	Y=X; % no change to input ringwidths
	loggy='xxxxxxxxxxxxxxx'; % 
	

end % of if k1
end % of while k1



