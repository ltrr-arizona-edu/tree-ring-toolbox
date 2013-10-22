function [u,yru,b,stats]=lwresid(xe,yrxe,xl,yrxl)
% dedicated function to get portion of latewood variance not predictable from earlywood

%--- Pull common period
yrmin = max([min(yrxe) min(yrxl)]);
yrmax = min([max(yrxe) max(yrxl)]);
L1 = yrxe>=yrmin & yrxe<=yrmax;
L2 = yrxl>=yrmin & yrxl<=yrmax;
ew = xe(L1);
lw = xl(L2);
mnlw = mean(lw);

if any(isnan(ew)) | any(isnan(lw));
   error('a NaN found in ew or lw');
end;

%--- Set year vector for result
yru =(yrmin:yrmax)';

%--- Build ones into predictor matrix
nsize  =length(yru);
X=[ones(nsize,1) ew];

%--- Do theregression
[b,BINT,R,RINT,stats] = REGRESS(lw,X,.05);
u = R + mnlw; % residuals, with mean set to mean of lw
