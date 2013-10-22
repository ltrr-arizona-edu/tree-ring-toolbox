function [yrs,yrf,sts,stf,chnd] = skprcs(rnb,st,m,flgop)
%
% USAGE : [yrs,yrf,sts,stf,chnd] = skprcs(rnb,st,m,flgop)
% A processing function for SKEMOD.M. Divides the SKE data 
% into 3 different blocks and returns only the first and 
% last blocks.
%
% INPUTS 
% ------
% rnb (1 x 1)	The input year index to be modified
% st (m x 2)	The 2-column old data matrix
% m (1 x 1)	Row size of the data matrix
% flgop (1 x 1)	A flag indicating if the user wants one
%		of the following : Add a LA ring, Delete a
%		a LA ring, Real to false, False to real
% OUTPUTS
% -------
% yrs (? x 1)	First block of year index
% sts (? x 1)	First block of SKE data
% yrf (? x 1)	Last block of year index
% stf (? x 1)	Last block of SKE data
% chnd (1 x 1)	Index corresponding to rnb
%
% NO USER WRITTEN FUNCTIONS NEEDED
%________________________________________________________

% Find the index corresponding to rnb
chnd=find(st(:,1)==rnb);
rnm=20;   	% Maximum deviation from the range of st
if isempty(chnd),
  while (st(1,1)-rnb)>rnm | (rnb-st(m,1))>rnm,
    rnb=input('Ring # out of range. Correct ring #  = ');
  end
  chn=find(st(:,1)<=rnb);
  if isempty(chn),
    chnd=1;
  else
    lchn=length(chn);
    chnd=chn(lchn);
  end
end

% Divide the data into three blocks
ls=st(:,1)<rnb;
sts=st(ls,2);
yrs=st(ls,1);
if flgop==2,
  lf=st(:,1)>(rnb+1);
  stf=st(lf,2); 
  yrf=st(lf,1)-1;
elseif flgop==1,
  lf=st(:,1)>=(rnb+1);
  stf=st(lf,2); 
  yrf=st(lf,1)+1;
elseif flgop==3,
  lf=st(:,1)>=(rnb+2);
  stf=st(lf,2); 
  yrf=st(lf,1)+1;
elseif flgop==4,
  lf=st(:,1)>=(rnb+3);
  stf=st(lf,2); 
  yrf=st(lf,1)-1;
end


% End of file
