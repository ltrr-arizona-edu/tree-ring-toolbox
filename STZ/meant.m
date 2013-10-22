function [IT,ITn,Tnms,Tyrs] = meant(IX,nms,yrs)
%
% DESCRIPTION : [IT,ITn,Tnms,Tyrs] = meant(IX,nms,yrs)
% Converts standard core indices, as culled by the mask from 
% BLDMASK.M, into tree indices by averaging all available core 
% indices for each tree in each year.
%
% INPUTS  :  IX (? x 1)   - Strung out vector of core indices
%	     nms (ns x 6) - Names for the cores
%	     yrs (ns x 3) - Years and starting row subscripts 
%			    for cores within IX
%
% OUTPUTS :  IT (? x 1)   - Strung out vector of tree indices
% 	     ITn (? x 1)  - Contain averaging info. Same length
%			    as IT.
%	     Tnms (? x 6) - Tree names
%	     Tyrs (? x 3) - Years and starting row indices for 
%			    trees in IT 
%______________________________________________________________



% Get parameters from TREENO.M.  This information does not yet take into
% account cores to be masked out though cmask.
[ti,Tnms,tn]=treeno(nms);
[ns,ms]=size(nms);
ntr=zeros(tn,1);

% Determine the number of cores in each tree
k=1; % initialize as 1 core for a tree
temn=1; % initialize as 1 tree in the data set
for i=2:ns,
  if ti(i,2)==ti((i-1),2); % same tree, another core
    temn=temn+1;
  else; % new tree, initialize as 1 core
    k=k+1;
    temn=1;
  end
  ntr(k)=temn;
end

yrbe=zeros(tn,2);
k=1;
for i=1:tn,
  yrbe(i,1)=min(min(yrs(k:k+ntr(i)-1,1)):...
             max(yrs(k:k+ntr(i)-1,1)));
  yrbe(i,2)=max(min(yrs(k:k+ntr(i)-1,2)):...
             max(yrs(k:k+ntr(i)-1,2)));
  k=k+ntr(i);
end

% Put yrbe matrix in Tyrs
Tyrs=[yrbe yrbe(:,2)-yrbe(:,1)+1];
   
% Initialize IT and ITn
itlen=sum(yrbe(:,2)-yrbe(:,1)+ones(tn,1));
IT=ones(itlen,1)*NaN;
ITn=ones(itlen,1)*NaN;

% Calculate the starting and ending row indices of IT in which the 
% tree indices for each tree will be stored
itin=zeros(tn,1);
sm=0;
for i=1:tn,
  sm=sm+yrbe(i,2)-yrbe(i,1)+1;
  itin(i)=sm;
end

% Loop over the tn trees to compute tree indices and put in proper 
% slots in IT
init=0;
k=1;
for i=1:tn,
  if ntr(i)==1,		% If tree has only one core
    IT(init+1:itin(i))=IX(yrs(k,3):yrs(k,3)+yrs(k,2)-yrs(k,1));
    tempit=ones(length(IT(init+1:itin(i))),1);
    % Put NaN's in ITn vector if IT has NaN
    tempit(isnan(IT(init+1:itin(i))))=ones(sum(isnan(IT(init+1:...
		itin(i)))),1)*NaN;
    ITn(init+1:itin(i))=tempit;
  else
    % Initialize the matrix of core indices for a single tree
    temit=zeros(yrbe(i,2)-yrbe(i,1)+1,ntr(i));
    % Fill the temit matrix with appropriate core indices from IX
    for j=1:ntr(i),
      temit(yrs(k+j-1,1)-yrbe(i,1)+1:yrs(k+j-1,2)-yrbe(i,1)+1,j)=...
         IX(yrs(k+j-1,3):yrs(k+j-1,3)+yrs(k+j-1,2)-yrs(k+j-1,1));
    end
    % Build the ITn vector from temit
    for j=1:yrbe(i,2)-yrbe(i,1)+1,
      lgcit=temit(j,:)~=0;
      ITn(init+j)=sum(lgcit');
    end
    % Fill the IT vector with averaged tree indices
    IT(init+1:itin(i))=sum(temit')'./ITn(init+1:itin(i));
  end 
  init=itin(i);
  k=k+ntr(i);
end
    

% End of file
