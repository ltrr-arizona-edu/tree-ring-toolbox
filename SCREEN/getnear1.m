function [iw,w,d,ny]=getnear1(NY,W,npref,D,dcrit)
% [iw,w,d,ny]=getnear1(NY,W,npref,D,dcrit)
% Get the column index of the climate station best satisfying
% criterion of nearness and number of years for modeling
% with the tree-ring series
%
%***************** IN ARGS **********************
%
% NY (m1 x n1)i number of years of overlap in previous output-
%		error modeling of chronologies with n1 nearest climate
%		stations.  Total of m1 chronologies.  NY from run of 
%		screen.m
% W (mW x nW)i  column index to nW nearest stations to each
%		chronology
% npref (1 x 1)i preferred minimum number of years for 
%		output-error modeling
% D (mD x nD)r distance (km) corresponding to W
% dcrit (1 x 1)r  critical distance which along with NY entry
%		allows station to be be classified "first class" for
%		output error modeling (see notes)
%
%
%******************* OUT ARGS ************************
%
% iw (mi x 1)i class of chron/climate pair (see below)
%		1=first class
%		2=second class
% w (mw x 1)i  selected climate station for each chronology
%
%
%****************** NOTES *******************************
%
% getnear1.m written to go along with screen.m, the screening of
% chrons vs climatic series.  Previous run of screen.m models
% is needed to get NY.  W and D also would have been computed
% ealier (see screen.m).  The initial run of screen.m modeled
% each chron against each of the nearest 5 (say) stations.  
% That run also indicated how much overlap (years) between each
% chron and seasonal climate series -- info now in NY.  Now we
% want to reduce the problem and repeat the modeling, but only
% for one station climate series paired with each chron.  The
% pairing is done so that there are two types of pairs:
%
% First class:  nearest climate station satisfying (1) at least
%		npref years of overlap data for modeling and (2) dcrit
%		km or closer to the tree site
% Second Class:  No stations in W satisfy the two above 
%		criteria.  Resort to the nearest station with an 
%		"acceptable" data overlap.  "Acceptable" is computed as
%		min(max(NY')), so that every chron will be modelable.
%
% The output column vector w should be stored in the appropriate
% screen.m input .mat file.  For example, in sip1104a.mat, after
% running getnear1.m


a = NaN;
a1=1;
a2=2;
[m1,n1]=size(NY);
iw=a(ones(m1,1),:);
w=a(ones(m1,1),:);
d=a(ones(m1,1),:);
ny=a(ones(m1,1),:);


% compute threshold number of required years.  
nthresh=min(max(NY'));

% Compute row index to chrons that have satisfy npref years
L1= (NY>=npref & D < dcrit);
I1 =  find(any(L1')'); % row pointer to subset of npref chrons
s1=length(I1);


L2 = (NY>=nthresh);
I2 =  find(~any(L1')'); % row pointer to subset that will
%		have to resort to nthresh threshold

if ~isempty(I2);
	s2=length(I2);
	iw(I2)=a2(ones(s2,1),:);
end

if ~isempty(I1);
	s1=length(I1);
	iw(I1)=a1(ones(s1,1),:);
end





%**********************
% Work on the first-class chronologies
LF = L1(I1,:);
[m2,n2]=size(LF);
LF=[zeros(m2,1) LF];
diff1 = (diff(LF'))';

% From right to left so that nearest is final choice

w1 = a(ones(m2,1),:);
d1 = a(ones(m2,1),:);
ny1 = a(ones(m2,1),:);
W2=W(I1,:);
D2=D(I1,:);
NY2=NY(I1,:);

for n = n1:-1:1;
	j = diff1(:,n);
	jj = j==1;
	jf=find(jj);
	js = sum(jj);
	if js>0;
		w1(jj)=W2(jf,n);
		d1(jj)=D2(jf,n);
		ny1(jj)=NY2(jf,n);
	end
end
w(I1) = w1;
d(I1) = d1;
ny(I1) = ny1;


%****************************
% NOW ON SECOND CLASS
% Accept nearest station with at least nthresh years for modeling

LF = L2(I2,:);
[m2,n2]=size(LF);
LF=[zeros(m2,1) LF];
diff1 = (diff(LF'))';

% From right to left so that nearest is final choice

w1 = a(ones(m2,1),:);
d1 = a(ones(m2,1),:);
ny1 = a(ones(m2,1),:);

W2=W(I2,:);
D2=D(I2,:);
NY2=NY(I2,:);

for n = n1:-1:1;
	j = diff1(:,n);
	jj = j==1;
	jf=find(jj);
	js = sum(jj);
	if js>0;
		w1(jj)=W2(jf,n);
		d1(jj)=D2(jf,n);
		ny1(jj)=NY2(jf,n);
	end
end
w(I2) = w1;
d(I2) = d1;
ny(I2) = ny1;
