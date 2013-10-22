function [rim,rsm,lnt,lny]=intprs1(r,YM,yrY,lxs)
%
% USAGE : [rim,rsm,lnt,lny]=intprs1(r,YM,yrY,lxs)
%   Processes the data passed from PROCES1.M
%
%
% INPUTS
%-------
% r (ny x 1)	A vector containing statistical info 
%		from Pearson/Spearman, Sign, or 
%		Hypergeometric tests
% YM (my x ny)	Time series matrix built out of vector y
%		by shifting its index by one
% yrY (ny x 1)	Year vector corresponding to end 
%		elements of each column of YM
% lxs (my x 1)	Length of X segment
%
%
% OUTPUTS
%--------
% rim (3 x 1)	Index vector to r having highest 3 values
% rsm (3 x 1)	Vector containing 3 highest values of r
% lny (my x 3)	Matrix containing 3 column vectors in YM
%		corresponding to 3 rim values
% lnt (my x 3)	Year matrix corresponding to lny
%
%
% NO USER WRITTEN FUNCTION NEEDED
%_____________________________________________________

% Sort r in ascending order
[rs,ri]=sort(r);
rln=length(rs);

% Pick out 3 top values from rs
rsm=rs(rln:-1:rln-2);
rim=yrY(ri(rln:-1:rln-2));

% Store the 3 YM columns and their corresponding 
% year vectors in matrices
lny=[];lnt=[];
for i=1:3,
  ny=YM(:,find(yrY==rim(i)));
  nt=(rim(i)-lxs+1:rim(i))';
  lny=[lny ny];
  lnt=[lnt nt];
end

% End of file
