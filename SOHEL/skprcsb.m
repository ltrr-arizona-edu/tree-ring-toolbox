function [yrs,yrf,sts,stf,chnd] = skprcs(rnb,yvo,st,m,flgop)
%
% USAGE : [yrs,yrf,sts,stf,chnd] = skprcs(rnb,yvo,st,m,flgop)
%
%	A processing function for SKEMOD.M
%_____________________________________________________________________

chnd=find(yvo==rnb);
chd=chnd;
if isempty(chnd),
  while (yvo(1)-rnb)>20 | (rnb-yvo(m))>20,
    rnb=input('Ring # out of range. Correct ring #  = ');
  end
  chn=find(yvo<=rnb);
  if isempty(chn),
    chnd=1;
  else
    lchn=length(chn);
    chnd=chn(lchn);
  end
end

if rnb>yvo(m),
  sts=st;
  stf=[];
  yrs=yvo(1:m);
  yrf=rnb+1;
elseif rnb<yvo(1),
  sts=[];
  stf=st; 
  yrs==[];
  if flgop==1, 
    yrf=yvo(1:m)+1;
  elseif flgop==2,
    yrf=yvo(1:m,1)-1;
  elseif flgop==3,
    if (rnb+1)<yvo(1),
      yrf=yvo(1:m)+1;
    else
      yrf=yvo(2:m)+1;
    end
  end
else
  if isempty(chd),
    sts=st(1:chnd,:);
    stf=st(chnd+1:m,:);
    yrs=yvo(1:chnd,1);
    if flgop==2,
      yrf=yvo(chnd+1:m)-1;
    elseif flgop==1,
      yrf=yvo(chnd+1:m)+1;
    elseif flgop==3,
      if (rnb+1)<yvo(chnd+1),
        yrf=yvo(chnd+1:m)+1;
      else
        yrf=yvo(chnd+2:m)+1;
      end 
    end
  else
    sts=st(1:chnd-1,:);
    yrs=yvo(1:chnd-1,1);
    if flgop==2,
      if (rnb+1)<yvo(chnd+1),
        stf=st(chnd+1:m,:); 
        yrf=yvo(chnd+1:m,1)-1;
      else
        stf=st(chnd+2:m,:); 
        yrf=yvo(chnd+2:m,1)-1;
      end
    elseif flgop==1,
      stf=st(chnd+1:m,:); 
      yrf=yvo(chnd+1:m,1)+1;
    elseif flgop==3,
      lf=yvo>=(rnb+2);
      stf=st(lf,:); 
      yrf=yvo(lf)+1;
    elseif flgop==4,
      if (rnb+1)<yvo(chnd+1),
        stf=st(chnd+2:m,:); 
        yrf=yvo(chnd+2:m,1)-1;
      else
        stf=st(chnd+3:m,:); 
        yrf=yvo(chnd+3:m,1)-1;
      end
    end
  end
end


% End of file
