function [S,yrvn,xvn,xv1,nblk] = blocdat(hf1,scid,xind1,yrvn,xvn,yrv,xv,xv1,sgc,S)

   % Blocking out specified data segments
   nblk=menu('# Block seg. ?','1','2','3','4','5','Continue');

   %nblk=shlmnu('# Block seg. ?','1','2','3','4','5','Continue');
   if nblk==6, return; end
   nbch=menu('Input type ?','Graphical','Keyboard');

   if nbch==1,
    figure(hf1);
    jh1=jdisp('Click at the corner points of each block');
    pause(1);
    close(jh1);
    [bsgx,bsgy]=ginput(2*nblk);
    bsgx=round(bsgx);
   else
    jh1=jdisp('Please return to command window for input');
    pause(1);
    close(jh1);
    bsgx=zeros(2*nblk,1);
    for i=1:nblk,
      bsgx(2*i-1)=input(['Beginning year for block segment # ',int2str(i),' = ']);
      bsgx(2*i)=input(['Ending year block segment # ',int2str(i),' = ']);
    end
   end   

   xindb=zeros(2*nblk,1);
   for i=1:2*nblk,
    S(scid,12+i)=bsgx(i);
    if bsgx(i)<yrvn(1),
      bsgx(i)=yrvn(1);
    elseif bsgx(i)>yrvn(length(yrvn)),
      bsgx(i)=yrvn(length(yrvn));
    end
    xindb(i)=find(yrvn==bsgx(i));
   end

   % Fill blocked intervals of data vector and years vector with NaNs
   for i=1:nblk,
    lbdm=length(xindb(2*i-1):xindb(2*i));
    yrvn(xindb(2*i-1):xindb(2*i))=ones(lbdm,1)*NaN;
    xvn(xindb(2*i-1):xindb(2*i))=ones(lbdm,1)*NaN;
    xv1(xind1+xindb(2*i-1)-1:xind1+xindb(2*i)-1)=ones(lbdm,1)*NaN;
   end


% End of file
