function [I,T,n]=treenum(nms,cmask)
%
% Find out how many trees are in the data set,after masking out
%   unwanted cores
% Assign sequential tree number to each core; '0' for masked cores
% Make a string matrix of tree ids for all selected (non-masked) cores
%
% D Meko 6-22-95
%
%**********************  IN ARGS ******************
%
% nms -- core ids of all cores in the original .rwl file
% cmask -- col vector of 1 or 0 .  1 means include this core
%      in building trees;  0 means exclude this core
%
%
%********************  OUT ARGS *******************
%
% I (? X 2)i   tree-id numeric integer matrix. 
%    Each row holds info for a given core as ordered in nms
%    Cols defined as follows:
%	1-sequential number of core
%	2-sequential tree number (0 if a masked core)
%	
% T (n x ?)s  Unique tree ids for each of n trees to be included
%	in subsequent analysis;  masked cores do not contribute to
%	the tree count.  
% n (1 x 1)i  number of trees represented by the selected 
%	(non-masked) cores
%
%
%******************* USER-WRITTEN FUNCTIONS CALLED **********
%
% maskind.m -- uses the core-mask vector (cmask) to make a matrix
%	of indices for testing for finding a "new" tree in 
%	sequential stepping through the cores
%
%******************* METHOD **********************************
%
% User must follow a specific core-naming and core-numbering convention
% for treenum.m to work correctly.  See the general program
% instructions for the naming.
%
% The tree ids (T) are "re-built" from the core ids to ensure a 
% consistent tree-id format.  The 3-letter site code and numeric
% tree number is extracted from the core-id string vector.  The tree
% number is also extracted, and converted to a 3-char string.  The
% site code and tree-number are then concatenated into a tree id.
% The method handles inconsistencies in tree numbering in the input
% files.  For example, core names FED2A and FED02A do not yield tree
% ids FED2 AND FED02, but yield the single tree id FED2. 
%
% EXAMPLES OF SITE/TREE/CORE IDS ASSIGNED AS "TREE ID"
%
% pad1a --> pad1, site pad, tree 1
% pad01a --> pad1,  site pad, tree 1
% pad106a --> pad106, site pad, tree 106
% pad01a1 --> pad1, site pad, tree 1
% 407011 --> 4071,  site 407, tree 1
% 407012 --> 4071,  site 407, tree 1
% 40712  --> 4071,  site 407, tree 1
% 40722  --> 4072,  site 407, tree 2
%
% CAUTION.  Problem if id has no core identifier.  For example,
%  pad12  would be converted to pad1, which assumes the 2 is a
%  core id.  In fact, pad12 might be intended as pad tree 12.
%  I haven't yet run across this as a problem.

[ms,ns]=size(nms);
pref = nms(:,1:3);  % Site code assumed in cols 1-3

I=zeros(ms,3);
dd = logical(zeros(ms,1));

% Get the numeric part of id with the tree number
for i = 1:ms;
	txt = nms(i,:);
	txt=deblank(txt); % lop off trailing blanks
	txt(1:3) =[]; % lop off site code

	% Remaining part of name should have a numeric tree number
	% followed by either (1) a letter code for core, or (2) a
	% sequential numeric code for core

	% Is there a letter in the remaining part of name
	j = find(isletter(txt));

	% Warn if the 4th char in id is a letter
	if ~isempty(j) & isletter(j(1));  
		disp(['Core id ' nms(i,:)]);
		disp('Char 4 in above name a letter');
		disp('Press any key to continue');
		pause
	end
	
	% Only remaining letter should be the 1-char core identifier
	% Warn if more than 2 leters in part after site code
	if length(j)>1;
		disp(['Core id ' nms(i,:)]);
		disp('More than 2 letters after 3-char site code');
		disp('Press any key to continue');
		pause
		txt(j:length(txt)) = []; % lop off all remaining chars from
  		% the A, B, or whatever, to the end of string
	elseif length(j)==1;
		txt(j:length(txt)) = []; % lop off all remaining chars from
  		% the A, B, or whatever, to the end of string
	elseif isempty(j); % no letters in remaining part of name
		txt(length(txt))=[]; % strip off the last character, assumed
			% to be a numeric core sequence number
		if isempty(txt),
			error('Tree ID empty after stripping rightmost number');
		end
	else
		error('Unidentifiable core id format');
	end
	I(i,[1 3])=[i str2num(txt)];
end


% Find out how many different trees.  A change to a new tree is
% assumed when the tree number changes, or if the tree number 
% remains the same but the 3-letter code changes in sequential
% core ids.  Cores with cmask(i) set to 0 are skipped in this
% tabulation

% Get matrix of row indices for nonmasked cores and their first
% preceding non-masked core for comparison.  The number of rows in
% K is the number of non-masked cores minus 1.  The second col of K
% is the row index in nms of the second-on non-masked cores.  The
% first col of K is the corresponding preceding core for comparison.
% kgo holds the row index of the first nonmasked core in nms.
[K,kgo]=maskind(cmask);



B1 = I(K(:,1),3);
B2 = I(K(:,2),3);
Bd = B2 - B1; % elements will be zero if no change in tree number
Bd=Bd~=0;

C1 = pref(K(:,1),:);
C2 = pref(K(:,2),:);
[mc,nc]=size(C1);
for j = 1:mc;
   Cd(j)=strcmp(C1(j,:),C2(j,:));
   Cd=Cd;
end


d = Bd | (~Bd & ~Cd');


% Assign a sequential tree number to each core in nms

dd(K(:,2)) = d;
dd(kgo)=1;
n = sum(dd);

I(:,2)=cumsum(dd);


% Build string matrix of tree names
ttt=pref(dd,:); % 3-letter site codes
I2 = I(dd,3); % the tree numbers for all cores in nms
m = max(I2); % largest tree number
T1 = num2str(m); 
for i = 1:n;
   T1=str2mat(T1,num2str(I2(i)));
end
T1(1,:)=[]; % drop the initial row;

T = [ttt T1];

% Fill slots in I with zeros for masked cores
I (~cmask,[2 3])=zeros(sum(~cmask),2);
I(:,3)=[]; % I3 only useful for debugging
