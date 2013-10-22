function ms = meansen1(X,yrs)
% meansen1: mean sensitivity for multiple core or tree indices
% ms = meansen1(X,yrs);
% Last revised 9-2-99
%
% In tree-ring standardization via corei.m, computes mean sensitivity. 
%
% INPUTS  :  X (? x 1) - Strung out vector containing core indices
%  	     yrs (ns x 3) - Year matrix 
%
% OUTPUTS :  ms (ns x 1) - Mean sensitivity vector
%__________________________________________________________________

% Initialize ms
[ns,dms]=size(yrs);
ms=ones(ns,1)*NaN;

% Loop through for ms
for scid=1:ns,
 % Cull out individual core indices
 yrv=(yrs(scid,1):yrs(scid,2))';
 xv=X(yrs(scid,3):yrs(scid,3)+yrs(scid,2)-yrs(scid,1));

 % Check for segmentation
 lzn=isnan(xv);
 zn=xv;
 zn(lzn)=[];
 if ~isempty(zn),
  lzo=zeros(length(lzn),1);
  for k=1:length(lzn)-1,
    if lzn(k)~=lzn(k+1),
      lzo(k)=k;
    end
  end
  lz=lzo;
  lz(find(lzo==0))=[];

  if ~isnan(xv(1)),
    lz=[0;lz];
  end
  if ~isnan(xv(length(xv))),
    lz=[lz;length(xv)];
  end
   
  msi=0;
  ni=0;

  % Compute the mean sensitivities 
  for k=1:2:length(lz)-1,
   xvt=xv(lz(k)+1:lz(k+1));
   lxvt=length(xvt);
   %sum=0;
   if lxvt>=2,
    xvt1=xvt;
    xvt1(lxvt)=[];
    xvt2=xvt;
    xvt2(1)=[];
    %for i=2:lxvt,
    % sum=sum+abs((xvt(i)-xvt(i-1))/(xvt(i)+xvt(i-1)));
    %end
    msi=msi+sum(abs((xvt2-xvt1)./(xvt2+xvt1)))*lxvt*2/(lxvt-1);
    ni=ni+lxvt;
   end
  end 
  ms(scid)=msi/ni;
 end
end


% End of file
   

