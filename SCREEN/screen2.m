function [iw,w,d,ny,chow]=screen2
% Select 'best' climate station for detailed screen3.m modeling of tree-ring series 
% CALL:  [iw,w,d,ny]=screen2(NY,W,npref,D,dcrit);
%
%****** INPUT: EITHER PROMPTED OR FROM EXISTING .MAT FILES ********
%
% NY (m1 x n1)i number of years of overlap in previous output-
%		error modeling of chronologies with n1 nearest climate
%		stations.  Total of m1 chronologies.  NY from run of 
%		screen1.m
% W (mW x nW)i  column index to nW nearest stations to each
%		chronology; from previous run of screen0.m,
%		stored in dist?.mat
% npref (1 x 1)i preferred minimum number of years for 
%		output-error modeling; screen-prompted
% D (mD x nD)r distance (km) corresponding to W; from screen0.m 
%		output
% dcrit (1 x 1)r  critical distance which along with NY entry
%		allows station to be be classified "first class" for
%		output error modeling (see notes); screen-prompted
%
%
%*********** SAVED OUTPUT IN  DIST?.MAT **************
%
% npref (1 x 1)i preferred minimum number of years for 
%		output-error modeling
% dcrit (1 x 1)r  critical distance which along with NY entry
%		allows station to be be classified "first class" for
%		output error modeling (see notes); screen-prompted
% iw (mi x 1)i class of chron/climate pair (see below)
%		1=first class
%		2=second class
% w (mw x 1)i  selected climate station for each chronology
% d (md x 1)i  distance chron to station
% ny (mw x 1)i number of years in OE model period
% chow (1 x 1)i  criterion used for selecting best station 
%		Either 1 or 2 (see notes).  1 considers separation distance
%		and years of overlap, but not variance explained. 2 gives
%		priority to high variance explained, as long as some minimum
%		threshold of overlap is met.  Example might be a minimum
%		overlap of 36 yrs.  Note that the number of years overlap
%		needed for "class 1" designation might be higher -- say 50 yr.
%
%****************** NOTES *******************************
%
% Gets and saves the column index of the climate station best satisfying
% criteria.  Two alternatives: 
% (1) best station means combination of near the tree-ring site and with
%		long overlap of unbroken years with tree ring series.  
%		Chooses nearest station that has at least npref years of overlap.
% (2) best is that with the highest variance explained proportion from 
%		screen1.m, as long as the station has at least some minimum specified
%    number of years of unbroken overlap with the tree-ring series.  If
%		station with highest variance explained does not overlap enough, then
%		check the station with next highest variance explained for sufficient
%		overlap.  Eventually arrive at a station with sufficient overlap, even
%		if it has the lowest variance explained of the candidate stations from
%		screen1.m
%
% screen2.m was written to go along with screen1.m, the screening of
% chrons vs climatic series.  Previous run of screen1.m models
% is needed to get NY, the number of overlap years for each climate/tree pair.
% W and D also would have been computed
% ealier (see screen1.m).  The initial run of screen1.m modeled
% each chron against each of the nearest 5 (say) stations.  
% That run also indicated how much overlap (years) between each
% chron and seasonal climate series -- info now in NY.  Now we
% want to reduce the problem and repeat the modeling, but only
% for one station's climate series paired with each chron. T  The
% pairing is done according to one of the two criteria listed above.
%
% Screen2.m also classifies the selected climate/tree pair as first or second class.
% First class:   (1) at least npref years of overlap data for modeling and (2) dcrit
%		km or less separation distance
% Second Class: all other pairs 


% Help in file selection by specifying letter code for tree-ring set
char=input('Letter code for tree-ring matrix: ','s');
fn1=[char '*.mat'];


% Get NY from screen1.m output
txt1=[char 'P#OS1.MAT -- file with NY'];
[file1,path1]=uigetfile(fn1,txt1);
eval(['load ',path1,file1]);

% Get the distance info from screen0.m output
txt2=['dist' char '.mat -- file with distance info'];
[file2,path2]=uigetfile('dist*.mat',txt2);
eval(['load ',path2,file2]);

% prompt for npref and dcrit
clc
disp('You will now be asked for the criteria for determining')
disp('whether a site-station pair is first or second class')
disp('in terms of separation distance and number of years of')
disp('uninterrupted overlapping data')
npref=input('At least this many years good overlap: ');
dcrit=input('Separation distance (km) <= this: ');


% Prompt for how to choose the best climate station
clc
disp('You will now be asked for the criterion for choosing a best climate');
disp('station for modeling.  Method 1 says consider only the separation');
disp('distance and the years of overlap of climate and tree rings');
disp('Method 2 says take the station that gave the highest tree-ring ');
disp('variance explained in the previous modeling in screen1.m, as long');
disp('as some minimum threshold of overlap is met');
disp('Note that the minimum threshold is automatically calculated from the');
disp('input matrix NY so that each tree series will have at least one');
disp('qualifying climate station');
disp(' ');
chow = input('Method for selecting best station:  ');



a = NaN;
a1=1;
a2=2;
[m1,n1]=size(NY);
iw=a(ones(m1,1),:);
w=a(ones(m1,1),:);
d=a(ones(m1,1),:);
ny=a(ones(m1,1),:);

% Make a matrix same size of NY with repeated rows of sequential col numbers
% So, if n1==5, AA would be [1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5; etc ]
AA  = 1:n1;
AA = repmat(AA,m1,1);

c = (1:m1)';

% compute default threshold number of required years.  
ndef=min(max(NY')); % Ensures that at least one station qualifies per chronology

% Set threshold to default if nthresh does not already exist
if ~exist('nthresh');
   nthresh = ndef;
end

% Allow overriding of default number of required years.  For example, automatic
% threshold of 36 yr might rule out great nearby stations with 30 yr overlap
prompt={'Enter the minimum acceptable number of yr of overlap for modeling: '};
def={num2str(nthresh)};
%def = {nthresh};
title='Allow to override automatic threshold';
lineno=1;
answer=inputdlg(prompt,title,lineno,def);
nthresh=str2num(answer{1});
if nthresh>ndef;
   error('Must pick an override lower than default nthresh');
end
 
% Make logical pointer to class 1 stations and class 2 stations
L1 = NY>=npref & D<=dcrit; % class 1 stations
L2 = ~L1; % class 2 stations

% Fill matrix same size of NY with 1 or 2 depending on series class
CL = a(ones(m1,1),ones(n1,1));
CL(L1)=1;
CL(L2)=2;

% Make logical pointer to stations with minimum acceptable year coverage
L3 = NY>=nthresh;
% Check that each row of NY has at least one acceptable year coverage
if ~all(any(L3'))';
   error('All rows of NY must contain at least one station with yr coverage');
end

% Get a cv indicating for each chronology, the station (col of W, D,etc) 
% for subsequent detailed modeling
if chow == 1; % if selecting nearest station with at least nthresh yr overlap
   % Find nearest station with at least nthresh yr coverage
   BB=AA .* L3;
   BB(BB==0)=m1+1; % do this so next statement will ignore zeros in finding minimum
   i3 =   (min(BB'))'; % cv of length m1 indicating column of selected series
   %   in NY, D, W for each row of NY, D, W
else; % select the station with at least nthresh yr coverage that 
   %explained the
   %highest tree ring variance in screen1.m analysis
   BB = S  .* L3; % BB will have zeros for stations with <nthresh yr overlap
   % Sort rows of BB
   [Y,I]=sort(BB');
   I=I'; % I now same size as NY, etc; right most col has col index for 
   % highest variance explained station
   i3 = I(:,n1);
end

%Convert row subscript c and col subscript i3 to linear indices
igrab = sub2ind([m1 n1],c,i3);

% Pull the station info for selected key station
CLvect = CL(:);
Wvect=W(:);
Dvect=D(:);
NYvect=NY(:);
iw = CLvect(igrab); % class of pair (1 or 2)
w = Wvect(igrab); % index %to climate station in mother matrix
d=Dvect(igrab); % distance (km)
ny=NYvect(igrab); % number of yr of overlap of climate and tree-ring series
   

% store results
txt4=[char 'P#OS2.MAT --store screen2 results here'];
[file4,path4]=uiputfile('*.mat',txt4);
eval(['save ',path4,file4, ' npref dcrit iw w d ny chow nthresh'])
