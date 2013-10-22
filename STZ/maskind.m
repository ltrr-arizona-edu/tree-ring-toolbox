function [I,igo] = maskind(cmask)
%
% Utility function to prepare matrix of comparison indices for 
% treenum.m.
%
%

J = find(cmask==1); % row subscripts of non-masked cores in nms
igo = J(1);
j2 = J(2);  % second non-masked core
nj = length(J);  % number of non-masked cores
lag = zeros(nj-1,1);

JJ = J(2:nj); % subset of row subscripts, 2nd nonmasked core on
for j=1:length(JJ);  % loop over nonmasked cores, beginning with second
   k = JJ(j);
   nz = 0;  % initialize number of masked cores preceding this core
   for i = k-1:-1:1; % loop backwards from preceding core
	if cmask(i)==0;
		nz = nz+1;
	else
		break
	end
   end
   lag(j) = nz+1;  % subscript of comparison core -- first
	% preceding non-masked core to check that same or differet
	% tree.
end
I = [JJ-lag JJ];


