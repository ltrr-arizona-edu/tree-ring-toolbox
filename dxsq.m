% dxsq.m  index sequential method of superposed epoch analysis
%
%  Programmed by D. Meko,  11-20-91
%  Latest changes: 2-12-91
%
% Given a series of key dates and a time series, computes (1) the average
% value of the time series at specified lags from the key dates (2) 95%
% confidence limits for these averages, and (3) probabilities of 0, 1, 2,
% 3,...  individual lag means exceeding the confidence limits.  
%
% Modified 5-21-91 to give output data array data.dat for use in designing
% publicn quality plot in GRAPHER

% Results are given in plots and tables.
%
%**********************  PRELIMINARIES ******************************
%
%  Create a .mat array holding the following arrays
%
% 		K - row vector with key years
%		X -  array of time series.  First column must be year.
%			Pgm will use one of the remaining columns for the
%			analysis.  The scalar w will denote which.  
%		L is # of lags.  For example, L=2 means use  lags
%				-2, -1,  0,  +1,  +2
%		W denotes which column of X to use as time series.
%			   	Year column doesnt count.  So W=1 means use
%				column 2 of X.
%		labs - 7-row string array of titles, axes labels, and
%			gtext.  Generated with dxsqlabs.m  
%		
%
%  Load the .mat array before running this pgm.
%
%
%***********************  list of variables *************************
%
% X (m x n) 	Input data array, each column a time series
% 	      	First column is year.
% K (1 x n) 	Input vector of key years.
% K1 (length(K) x 2*L+1   index for lagged values of x relative
%         to key years.  Index shifted forward by number of years in X.
%
% L 	number of years before and after key year to compute probabilities
% W (1 x n)   which time series of X to treat (1=column 2)

% ********* obsolete code for automatic prompting for file name ***
%
% file1=input('NAME OF .MAT STORAGE FILE: ','s');
%eval(['load ',file1]);
% disp('here')
%
% **********  end of obsolete code

kgo=zeros(length(K));     % initilize arrays for start and end rows
kstop=zeros(length(K));	  % for averaging 

clear Xthis
[m,n]=size(X);
yr=X(:,1);
Xthis=X(:,W+1);

%*************** Plot time series with key years marked **********
clg
hold off
plot(yr,Xthis,'r-');
hold on
xstd=std(Xthis);
xbar=mean(Xthis);
xdiff=xbar-xstd;
level=xbar(ones(m,1),:);
plot(yr,level,'g-');
plot(K',xdiff(ones(length(K),1),:),'w+');

gtext(labs(1,:));   % text for title on time series plot
xlabel(labs(2,:));
ylabel(labs(3,:));
gtext(labs(4,:));


pause
hold off 


Xbig=[Xthis;Xthis;Xthis];   % stack three copies of time series back to back

xxm=zeros(length(yr),2*L+1);
K3=zeros(length(K)*m,2*L+1);

% Form array of relative indices for key years

KTEMP= K - yr(1)+1 + m;  % shifted by  number of years in time series
KTEMP=KTEMP'* ones(1,2*L+1);    % col vector of key year indices

MTEMP= [-L:L];
MTEMP= ones(length(K),1) * MTEMP ;
 
K1= MTEMP + KTEMP;

%  Compute means at various lags

d=Xbig(K1);
e=ones(length(K),2*L+1);
e(:)=d;
xm=mean(e);   % means at various lags


% Form array of indices for "synthetic" series

K2=K1;   % initialize K2 for loop
oneblk = ones(length(K),2*L+1);

ksize=length(yr) *  length(K);   % # rows in g
kstop = [length(K):length(K):ksize];  % stop row for each synthetic
kgo  = kstop - (length(K)-1);    % go row for each synthetic series

disp('Patience, results will come.');

for j = 1:m
	K3(kgo(j):kstop(j),:) = K2 + oneblk;
	K2=K2+oneblk;
end

% Values of synthetic series at lags

f=Xbig(K3);
g=ones((m)*length(K),2*L+1);
g(:)=f;

%***********  form series of averages over length(K) years *****

for i = 1:length(yr)    % loop thru each synthetic series
	xxm(i,:) = mean(g(kgo(i):kstop(i),:));  % means over each
      %  set of snythetic key years 
end

%***************** compute selected percentile cutoffs *************

xsort=sort(xxm);

f1=[.05 .95];
f2=length(yr) .* f1;  %  index of 5th and 95 th percentile
c1=floor(f2);  % nearest integer below
c2=ceil(f2);
d1= xsort(c1);  
d2= xsort(c2);

inc=((f2-c1)/(c2-c1))* (d2-d1);  %delta displacement from lower value
cut=inc+d1;  % time series value of cutoffs

%**************** plot means at various lags, with 95% error bars ***

cutrows= cut' * ones(1,2*L+1);   % expand cutoffs into rows 

lag=-L:L;

[x1,y1]=bar(lag,xm);


bsize=length(x1);
base=mean(xsort(:,1)) .* ones(1,bsize);
down=cutrows(1,1) .* ones(1,bsize);
up=cutrows(2,1) .* ones(1,bsize);


for jj=1:length(y1),
	if y1(jj) == 0  
		y1(jj)=base(1);
	end
end


clc, home;

TBL1=[lag;xm;cutrows(1,:);cutrows(2,:)]';
disp('   LAG    MEAN          95% LIMITS');
disp(' ');
disp(TBL1);
pause,clc, home;

plot(x1,y1,'-r',x1,base,'-g',x1,up,'--r',x1,down,'--r')

data=[x1 y1 base' up' down'];

save data.dat  data /ascii

xlabel('LAG')

% label 95% confidence above bottom bar, 2/3 of way to right

ltemp=length(x1);
ltemp=fix(2*ltemp/3);
text(x1(ltemp),down(1,1),'95 % CL')


gtext(labs(5,:));
ylabel(labs(6,:));
gtext(labs(7,:));

gtext(['n1 = ',num2str(length(yr))]);
gtext(['n2 = ',num2str(length(K))]);


pause

% Compute number of simulations with 1,2,... etc of the means greater
% than the upper and lower than the lower confidence limits


%*****************  Analysis based on Upper 95 % Confidence Limit *****

D= (xxm >= up(ones(m,1),1:2*L+1)); % test all entrirs of xxm for
				% exceedance of upper 95% threshold
E = sum (D');  % row vector, number of lag-means exceeding threshold
				% in each year.  Max possible of 2*L+1.
E1 = E(ones((2*L+2),1),:);  %  duplicate E into 2*L+2 rows
F = [0:2*L+1];  % could be min of zero and max of 2*L+1 exceedances
		% in a given year
F1 = F(ones(m,1),:)';  % dupe F into m rows
G1= F1==E1;
H1 = sum(G1');
TBL2=[F;H1]';

clc, home;
disp('NUMBER OF SIMULATIONS WITH 0,1... MEANS GREATER THAN UPPER LIMIT')
disp(' ')
disp(TBL2);
pause,clc,home;

% convert number of simulations outside limits to empirical
% probabilities

M = sum(H1);   % should be equiv to total number of simulations
M1 = M(ones(2*L+2,1),:)';
M2 = H1 ./ M1;     % convert freqs to relative frequencies
TBL3 =  [F;M2]';

disp('EMPIRICAL PROBABILITY OF 0,1,... MEANS GREATER THAN UPPER LIMIT')
disp(' ');
disp(TBL3)
pause,clc,home;

% compute prob of less than or equal to 0,1,2... exceedances

N = cumsum(H1);
N1 = N ./ M1;
TBL4 = [F;N1]';

disp('CUMULATIVE PROBABILITY OF 0,1,2... MEANS GREATER THAN UPPER LIMIT')
disp(' ')
disp(TBL4);
pause, clc, home;


%***************  Analysis based on Lower 95% Confidence Limit ************

% Compute number of simulations with 1,2,... etc of the means less
% than the lower confidence limits

D= (xxm <= down(ones(m,1),1:2*L+1)); 
				
E = sum (D');  % row vector, number of lag-means less than threshold
				% in each year.  Max possible of 2*L+1.
E1 = E(ones((2*L+2),1),:);  %  duplicate E into 2*L+2 rows

G1= F1==E1;
H1 = sum(G1');
TBL5=[F;H1]';

clc, home;
disp('NUMBER OF SIMULATIONS WITH 0,1... MEANS BELOW LOWER LIMIT')
disp(' ')
disp(TBL5);
pause,clc,home;

% convert number of simulations outside limits to empirical
% probabilities

M = sum(H1);   % should be equiv to total number of simulations
M1 = M(ones(2*L+2,1),:)';
M2 = H1 ./ M1;     % convert freqs to relative frequencies
TBL6 =  [F;M2]';

disp('EMPIRICAL PROBABILITY OF 0,1,... MEANS BELOW LOWER LIMIT')
disp(' ');
disp(TBL6)
pause,clc,home;

% compute prob of less than or equal to 0,1,2... exceedances

N = cumsum(H1);
N1 = N ./ M1;
TBL7 = [F;N1]';

disp('CUMULATIVE PROBABILITY OF 0,1,2... MEANS BELOW LOWER LIMIT')
disp(' ')
disp(TBL7)
