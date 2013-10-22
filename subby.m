iLT=find(LTRI);
LTRI1=LTRI;
LTRI1(iLT(1)-maxar:iLT(1)-1)= ones(maxar,1);
LTRI1(iLT(1)-maxar:iLT(1)-1) = ones (maxar,1);
ysub=TRI(LTRI1,iw);
[eyar,arorder,vrat,arcs]=whit1(y,maxar,1);
tree=eyar(maxar+1:length(eyar));

