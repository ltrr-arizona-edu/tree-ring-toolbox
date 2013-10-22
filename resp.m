% resp.m  

% January, 1992, by D. Meko

% Application:  Screening of a tree-ring index against monthly
% precipitation and temperature data to find out the best seasonal 
% grouping of months for a climate signal

% Eight output matrices are optionally produced.  The tree-ring index
% is optionally prewhitened for each analysis.  Each matrix contains
% product-moment correl coefs between pairs of series of different
% types of data:

% 1. 	tree vs ppt -- F
% 2.  tree vs temp -- G
% 3.  tree vs pptW -- U
% 4.  tree vs tempW -- V
% 5.  ppt vs temp -- H
% 6.  pptW vs tempW -- S
% 7.  tree.ppt vs temp.ppt -- R
% 8.  tree.pptW vs temp.pptW -- T
% 9.  first order autocorrel coefs for ppt -- P
% 10. first order acfs for temp -- Q


% In the above, W indicates prewhitened by ar model, and y.x means
% residual of linear regression of y against x.  Each array
% (12 x 12) summarizes correlation analysis between 144 pairs of series. 

% Bar charts giving means and standard deviations for montly ppt and
% temp are also produced.

% Each correlation array is 12 x 12.  A given row represents a 
% specific ending month of a monthly grouping.  A given column 
% represents a specific numbr of months in the grouping.  

% Arrays are optionally saved as ASCII files in form directly usable
% by surfer for getting high-quality contoured plots.

%*************   PRELIMINARY STEPS  ****************************************

% Run FIL999.F  and REAR1.M  to get the ppt and temperature data into MATLAB
%  array form.  Each array must have 13 numeric-only cols.  These arrays
%  must not have header lines.  First col is the year;  remaining cols are
%  the monthly values for Jan through Dec.

% Program doesn't handle missing values.  MATLAB also won't 
% recognize arrays with blank elements.  So put some dummy value in for missing
% values in arrays beforehand.  The dummy value must be numeric--like
% 09999.

% Program will prompt for period on which you want to run the correlation
% analysis, which can be subset of years of entire arrays A and B.  Make
% sure all the data in the period for analysis is valid (i.e., does not
% contain special missing-value codes, like "99999").

% The tree-ring array is assumed to be at least 2 columns.  First column
% holds the year, succeeding cols the tree-ring indices for one or more
% cores or sites.  As for the PPT and temperature arrays, no alphnumerics
% are allowed, no header lines, and no blank cells.

%*************** PRELOADS  ****************************************

% X.MAT	monthly precip array and related files in matlab form
%	A	the precip array, values scaled by factor of 0.01 in pgm
%	B	the temp array, values scaled by factor of 0.1 in pgm
%	C   the tree-ring array
%  nhi - highest ar-order to consider fitting to tree rings or climate
% 	endmo - ending month of growth year

%***********   RESTRICTIONS   ************************************

% A, B must be of same size and cover same years
% 
% C  	can cover any period overlapping A and B, but must have at least 
% 	nhi more years at the start to allow for loss of data when 
%	whitening.  
%
% C 	cannot have missing data in the period specified for analysis
% 	

%***************   OUTPUT    *******************************************

% F	mtx of correl for trees against PPT
% G   mtx of correl for trees against Temp
% H   mtx of correl for PPT vs Temp
% P   mtx of first order autocorr coefs for PPT
% Q   mtx of first order autocorr coefs for Temp
% R   mtx of partail correls (resids of tree vs ppt) vs (resids of temp vs ppt)
% U   mtx of correl for trees vs whitened ppt
% V   mtx of correl for trees vs whitened temp
% S  like H, except ppt and temp prewhitened
% T  like R, except ppt and temp prewhitened 

% Col i in the above arrays holds an "i month" sum or mean
% Row j in the above arrays hold a sum or mean ending in a particular
%   month relative to a specified ending month for the tree-ring year.
%   For example, if endmo=8 (August), row 1 is for all monthly groupings
%   ending in August of year t;  row 2 is for groupings ending in July
%   of year t ..., and row 12 is for groupings ending in September of 
%   year t-1.



%****************   SELECTED VARIABLES  (EXCEPT OUTPUT ARRAYS)   *********************

% A (m1 x 13) the precip array
% B (m1 x 13) the temp array
% C (m2 x n2)  the tree-ring array
%
% endmo - ending month of growth year (8=august)
% d	(scalar) which tree-ring series to analyze--in col d+1 of C
% F1 - F strung out as a column vector
% F2 - x,y,F1 stored for saving as ascii array
% k - miscellaneous screen input: 
%   - save ascii files for surfer?
% k1 - if loop on scaling A,B, creating D, E
% k2 - for loop for treating a given tree-ring series
% k3 - code for which type of array to generate.  This code number
%    is suffixed to ascii file name of file F?.dat.  So, if array type
%    is 
% L1 (lcv) marks approp rows of D,E for this tree-ring series analysis
% L2 (lcv) marks approp rows of C
% m1,n1  of ppt array A and temperature array B
% m2,n2  of tree-ring array C
% n3     of series designator
% vect1 (1x12) ending month vector to plot along left side of output
%  	corr matrices.  Example for endmo=8: 
%	vect1=[9 10 11 12 1 2 3 4 5 6 7 8]
%
% x (cv) x-coords for elements of F, strung out as single column
%	Used by surfer.
% y (cv) y-coords (see x)


%********  LOADS, PRELOADS, AND PREALLOCATES   *****************

nhi=5;  % Temporarily hard-coded max ar model to use
kind=['F' 'G' 'H'  'P' 'Q'  'R' 'U' 'V' 'S' 'T']';
ty = ['D' 'E' 'D' 'E' 'D' 'E'];
tit=['TREE    VS PPT                  ' 
'TREE    VS TEMP                 '
'AR(1)     PPT                   '
'AR(1)    TEMP                   '
'PPT     VS TEMP                 '
'TREE.PPT VS TEMP.PPT            '
'TREE VS WHITENED PPT            '
'TREE VS WHITENED TEMP           '
'PREWHITENED PPT VS TEMP         '
'PREWHITENED TREE.PPT VS TEMP.PPT'];

clear C1 F G H P Q R U V S T;

% Check that necessary input exists

t = [exist('A') exist('B')  exist('C')   ...
    exist('endmo') exist('nhi') ];

if t~= [1 1 1  1  1]
	disp('YOU FORGOT TO LOAD APPROPRIATE ?.MAT')
	keyboard;
end

[m1,n1]=size(A);
[m2,n2]=size(B);
n3=1;             % Number of series to analyze, change later to while lp
[m4,n4]=size(C);

					
% Preallocate six key output arrays

F=zeros(12); G=zeros(12); H=zeros(12); P=zeros(12); Q=zeros(12);
R=zeros(12);  U=zeros(12); V=zeros(12); S=zeros(12);  T=zeros(12);  
OP=zeros(12,12);  OT=zeros(12,12);  VP=ones(12,12);
VT=ones(12,12);



for k2=1:n3;  % loop for each tree ring series to be analyzed

	d=input('WHICH TREE-RING SERIES? ');
	a(1)=input('FIRST YEAR FOR ANALYSIS PERIOD: ');
	a(2)=input('LAST YEAR FOR ANALYSIS PERIOD: ');
	
	arord=0;    % initialize ar order for tree rings to zero


%  Check that climate arrays exist and have been scales.
%  Form the appropriate 3-year arrays.  Scale ppt and temp.
	
	if ~exist('D'); % A,B need to be scaled, D,E created
		A(:,2:13)=A(:,2:13) * 0.01;  % scale ppt
		B(:,2:13)=B(:,2:13) * 0.1;   % scale temp

		% build D, E

		D= [A(3:m1,1)  A(1:m1-2,2:13)  A(2:m1-1,2:13)  A(3:m1,2:13)];
		E= [B(3:m1,1)  B(1:m1-2,2:13)  B(2:m1-1,2:13)  B(3:m1,2:13)];

	end;  % if for scaling ppt, temp and creating lagged arrays

	% Calc desired rows of D, E, C for current tree-ring series

	L1 = (D(:,1) >= a(k2,1) & D(:,1) <= a(k2,2)); % approp years 
	L2 = (C(:,1) >= a(k2,1) & C(:,1) <= a(k2,2)); % approp years 


%****************   CHECK THAT YEAR RANGE IS CONSISTENT   *****************

if (C(L2,1) ~= D(L1,1));
	error ('L2 AND L1 SPECIFY DIFFERENT YEARS IN C THAN IN D')
end

check=C(L2,1);
d1=(check(1):check(length(check)))';
if d1 ~= check 
	error ('L2 and L1 DO NOT GIVE UNBROKEN STRETCH OF YEARS')
end


%  Pull off the correct tree-ring series and its correct rows for the
%  analysis.  Prewhiten it with an up-to-nhi order ar model, and store
%  the original and whitened versions.

	C1=C(L2,d+1);  % grab the correct tree-ring series	

% Model as ar(up-to-order-nhi) to get ar order.  This model will be
% built on years pointed to by L2.  The first arorder residuals will
% not be valid.  

	[C2,arord,varrat,arcs]=whit1(C1,nhi,1);  % Find order for ar model.

	if arord ~= 0;  %  Fit same order ar model, but extend time series

%	arord years on front end so that get valid resids for
%	full length of C2

		nhiset=arord;
		j=find(L2==1);
		L3=L2;
		L3(j(1)-arord:j(1)-1) = ones(arord,1);
		[C2,arord,varrat,arcs]=whit1(C(L3,d+1),nhiset,2);
		C2(1:arord)=[];	

	
	else;  % ar order is zero -- null model

%	Whit1 will have returned arord=0, C2= C1, and the
%	(un-used) ar coefs and std devs  in arcs.
		varrat = 1.0;
	end

 	disp(arord);
 	disp(arcs);
 	disp(varrat);
 	plot(C(L2,1),C1,C(L2,1),C2);
 	title('ORIGINAL AND WHITENED TRI');
 	pause


% C1 now holds the original tree-ring index.
% C2 holds the prewhitened index.

	k3=1;

	while k3~=12 

	k3=menu('TYPE OF ARRAY WANTED','TREE VS PPT','TREE VS TEMP',...
  	'AR COEF FOR PPT','AR COEF FOR TEMP', 'PPT VS TEMP',...
	'TREE.PPT VS TEMP.PPT',...
   'TREE VS WHITENED PPT',...
   'TREE VS WHITENED TEMP',...
   'PPT VS TEMP, PREWHITENED',...
   'TREE.PPT VS TEMP.PPT, PREWHITENED',...
   'SINGLE MONTH MEANS AND STD DEVS','QUIT')

	if k3==12, break, end;

%*************************************************************************


	if k3>=1 & k3 <=4;  %  want tree vs ppt or temp, or acs of ppt,temp
		eval([kind(k3),'=corr1(C1,',ty(k3),',L1,endmo,k3);']);

	elseif (k3 >= 5 & k3 <=10) ;  % want ppt vs temp or tree.ppt vs temp.ppt
%         possibly with prewhitening of climate
		kw=1;  % prewhiten, unless  change below
		C3=C2;  % use prewhitened tree 

		if k3==5 | k3==6 
			kw=0;
			C3=C1;
		end;

		eval(['[',kind(k3),',OP,OT,VP,VT]=corr3new(C3,D,E,L1,endmo,arord,k3,kw,nhi);']);

	elseif k3==11
		[xM,xS]=corr2(D,E,L1,endmo,k3);
		xaxe=endmo-11:endmo;
		subplot(221)
		bar(xaxe,xM(1:12));
		title ('MEAN TEMP');
		
		subplot(222)
		bar(xaxe,xM(13:24));
		title('MEAN PPT');
		
		subplot(223)
		bar(xaxe,xS(1:12));
		title('STND DEV TEMP');
		
		subplot(224)
		bar (xaxe,xS(13:24));
		title('STND DEV PPT');
		pause
	end
	
	if k3~=11
		lowr=eval(['min(min(',kind(k3),'))']);
		hir=eval(['max(max(',kind(k3),'))']);
	 
%		up=eval(['0.1 * floor(max(max(',kind(k3),'*10)))']);  % upper 
%		down=eval(['0.1 * ceil(min(min(',kind(k3),'*10)))']);  % lower 

		eval(['contour(',kind(k3),',.8:-.1:-.8,1:12,endmo-11:endmo)']);

		xlabel('NUMBER OF MONTHS')
		ylabel('ENDING MONTH ')
		title([tit(k3,:)]);
		text(.5,.25,['PERIOD ',int2str(a(1)),'-',int2str(a(2))],'sc');
		text(.5,.21,['MINIMUM R = ',num2str(lowr)],'sc');
		text(.5,.17,['MAXIMUM R = ',num2str(hir)],'sc');
		pause
	end  % if loop for k3~=11

%************  SURFER OUTPUT FOR FINE PLOTS OF CORR MTXS  ***********

k=input('MAKE SURFER FILES? Y/N [N]','s');
if isempty(k)
	k=='N';
end

if k=='Y'
	x=(1:12);
	x=x(ones(12,1),:);
	x=x(:);

	y=(endmo:-1:endmo-11)';
	y=y(:,ones(1,12));
	y = y(:);

	if k3==1
		F1=F(:);
	elseif k3==2
		F1=G(:);
	elseif k3==3
		F1=H(:);
	elseif k3==5
		F1=P(:);
	elseif k3==6
		F1=R(:);
	elseif k3==7
		F1=U(:);
	elseif k3==8
		F1=V(:);
	elseif k3==9
		F1=S(:);
	elseif k3==10
		F1=T(:);
	end
		F2=[x y F1];

	eval(['save c:\work\',kind(k3),'.dat F2 -ascii']);

end;  % k if for making ascii .dat file


end;   % k3  while for type of analysis

end;  % k2 for loop for each series 
