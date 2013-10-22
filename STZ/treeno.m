function [Inms,tnms,n]=treeno(nms)

% DESCRIPTION : [Inms,tnms,n]=treeno(nms);
% This function finds the number of tree in nms.
%
% INPUTS :  nms (ns x 6) -  String matrix containing the core
%             indices of the trees
%
% OUTPUTS : Inms (? x 2) - Contains the tree sequence # and the
%			   tree numbers in nms.
%	    tnms (? x 6) - A string matrix containing the 
%	                   tree names.
% 	    n (1 x 1)    - Total number of trees.
%____________________________________________________________

% Get the dimension of nms
[ms,ns]=size(nms);

% Initialize Inms and n
Inms=zeros(ms,2);
Inms(:,1)=(1:ms)';
n=1;

% Check for column dimension of nms
if isempty(str2num(nms(1,4:5))),
  ns=5;
else
  ns=6;
end

% Initialize tnms
tnms=nms(1,1:ns-1);

% Find the core index numbers
for i=1:ms,
  Inms(i,2)=str2num(nms(i,4:ns-1));
  if i>1, 
   if Inms(i,2)~=Inms(i-1,2),
    n=n+1;
    tnms=[tnms;nms(i,1:ns-1)];
   end
  end
end

 
% End of file