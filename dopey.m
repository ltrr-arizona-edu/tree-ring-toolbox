function datout=respfn1b(datin);
% respfn1b: subfunction of respfun1. PCR analysis and graphics for response function
% CALL: datout=respfn1b(datin);
%
% Meko 10-6-98
%
%******* IN
%
% datin{}
%  {1}=y; % standardized tree-ring series
%  {2}=Z; % standardized climate data(nseasv columns)
%  {3}=nseasv; % number of seasonsal variables = num of seasons times 2
%  {4}=svlab; % names of seasonalized variables
%  {5}=strpd; % string for period of analysis
%  {6}=treenm; % char name of tree-ring chron
%  {7} (1 x 2)i  options
%     1=kconf  method for confidence bands around PCR weights
%         1==statistical (t-dist)
%         2==bootstrap
%         3==circular shifted
%         4==monte carlo
%  {8}=fwind  first figure window available to this function
%  {9} Clr (3 x 9)r color matrix (see respfun1.m)
%  {10} dlab (2 x 1)s data-type labels (e.g., ['P';'T']
%  {11} slab (1 x ?)c  cell vector of season names
%
%
%********* OUT
%
% datout{}
%  {1} V eigenvectors of the climate variables, each col an eigenvector
%  {2} teach (? x 1)r percentage of climate variance for each PC
%  {3} V1a (? x 1)r decimal fraction of tree variance described by each
%      climate PC
%
%
%**********  NOTES
%
%
% Effective sample size is not adjusted for autocorrelation in the tree-ring
% or climate variables.  Keep that in mind when comparing the significance
% of weights on monthly climate by this method and by the correlation method
% in respfun1.  The correlation method does optionally adjust the effective
% sample size in computation of signifance.
%
%
% Three figure windows are produced.
% First window:
%  top plot: PCR weights in bar chart, with 99% and 95% signif color marked.
%   Ordering along x axis is by number of climate PC.  
% Second Window: R-sqd and RE plot vs step, with RE based on cross-validation
%     x-axis ordered as climate PC#. Bars annotated with entrystep 
% Third Window:
%  Response function weights on monthly climate variables, from PCR model selected
%     interactively by user after he looks at first window. Third window not produced
%     until user clicks on PC in plot in second window
%
% User actions.  User looks at plots in second window and decides a cutoff number
% of climate PCs to include as predictors in the model. User allowed only
% to use PCs whose regression coeffs signif at .05 or better. User can use fewer
% PC's. In limit, might just use first PC to enter.  User tells function the last
% PC he wants to include by pointing to it's bar on the top plot in window 1.



%------- UNLOAD
y = datin{1}; 
Z = datin{2};
nseasv = datin{3};
svlab=datin{4};
strpd=datin{5};
treenm=datin{6};
kconf=datin{7}(1);
kcolor=datin{7}(2);
fwind = datin{8};
Clr=datin{9};
dlab=datin{10};
slab=datin{11};

%----- Size
[mz,nz]=size(Z);

nullfg2=0;  % this will change to 1 if no significant predictors
[RR,V,L,S,F]=eigen1(Z,1);  % pca on correl mtx of Z
W=Z*V;  %  Amplitudes of PCs of Z (cols).  Recall that Z contains standardized
% variables, and so they have mean zero.  As a result, means of cols of W
% are also zero

% Compute variance explained of climate by each climate PC

ts=sum(diag(L));
teach=100 * diag(L) ./ ts;


a=W\y;  % Estim coefs of model to predict tree from clim eig amplitudes
ep = y - W*a;  % Compute regression residuals

% Compute pct variance of tree growth explained by each climate 
% eigenvector.

Rtemp=corrcoef([y W]);
V1a=(Rtemp(1,2:nz+1)) .^2;  % proportion variance expld
V1b= sum(V1a);  % Total proportion variance expld by model with all pc's

epsq=ep' * ep;  % sum of squares of residuals of regr of tree variable on the
			  % full set of  pc amps of climate array.

term1=epsq(ones(nz,1),:); % Terms needed for eq 8.85, p. 245 in 
aa=(mz-nz-1);             % Mardia et al. 1979
term2=aa(ones(nz,1),:);
term3= sqrt(mz * diag(L));

% Find .975 and .995 probability points of "t" distribution:  for 95%
% and 99% confidence limits around regression coefficients a.
% Only the 95% cl is used in the plots, but the 99% is there if you need
% it.
t1=tinv(.975,aa);
t2=tinv(.995,aa);
t95=((t1 * term1 ./ term2) .^2) ./  term3;
t99=((t2 * term1 ./ term2) .^2) ./  term3;

% Logical pointer to PCR stdzd regression weights signif at 99%, 95%, and not at 95%
Lc = logical(zeros(nseasv,3));
Lc(:,3) = abs(a) < abs(t95);
Lc(:,2) = abs(a) >= abs(t95) & abs(a)<t99;
Lc(:,1) = abs(a)>=t99;



J = abs(a) < abs(t95);  % Ones point to elements to zero out
Jr = ~J;  % 0/1 cv, 1s point to nonzero elements of a
ev = epsq / mz;  %  Variance of residuals from PC regression using all
%		PCs as predictors.
vexp=sum(V1a(Jr));  % Total proport. variance tree growth explained by the 
%	restricted set of climate eigenvectors
Jrf = find(Jr);  % Subscripts of cols of A corresp to restricted set
%      of PC predictors

%******  F-level and significance for PC regression using only
%        coefs sig at 95% level



if sum(Jr)==0, nullfg2=1; end;  % No significant PC coefficients--null model

% F-test for significance of the multiple correlation coefficient
% See Panofsky and Brier, p. 113.  In table 32 of that reference, see that
% the F is computed from R-squared, N and p.  R-squared is the proportion of
% variance accounted for by regression, which in my code is vexp.  N is the number
% of observations, which is mz; and p is the number of predictors, which is 
% sum(Jr) in my notation.


Lcrit = logical(zeros(3,nseasv));

if nullfg2==0;    %  At least one signif pc coefficient; proceed
   dfreg=sum(Jr); % degrees of freedom for regression sum of squares
   dfres=mz-sum(Jr)-1;  % deg freedom for residual sum of squares
   F2=vexp*dfres/((1-vexp)*dfreg); % Compute F ratio
   ff2 = finv(0.95,dfreg,dfres); % table value of F
   
      
   %***********  TRANSFORM THE PC-REGRESSION WEIGHTS BACK INTO 
   %***********  WEIGHTS ON THE INDIVIDUAL MONTHLY CLIMATE VARIABLES, AND
   %***********  COMPUTE ASSOCIATED ERROR BARS (MARDIA ET AL. 1979, P. 246)
   
   % Retain only those components whose  tree vs pc  regression coefs were
   % significant at 95% level.  Form a restricted least squares estimator
   % by setting all non-significant a's to zero.
   
   arest=a;  % initialize restricted estimator
   arest(J)=zeros(sum(J),1);
   b=V * arest;  %   Transformed coefs -- on the original monthly climate
   
   
   bv=zeros(length(Jr),1);
   %bv=zeros(Jr);  % Preallocate
   
   Ld=diag(L);
   
   for i=1:nz;  % for each monthly climate variable
      bv(i) = (ev / mz) * sum ((V(i,Jr) .^2) ./ (Ld(Jr))');
   end
   
   bvse=sqrt(bv);  % standard errors for coefs of transformed
   %		model.
   
   % Call matlab stat function norminv to get level of coefs different from
   % zero at the 99% and 95% levels Mardia, p. 246
   crit99=norminv([.995 ],0,bvse');
   crit95=norminv(.975,0,bvse');
   % Logical pointer to correl coefs signif at 99%, 95%, and not at 95%
   Lcrit(3,:) = abs(b') < abs(crit95);
   Lcrit(2,:) = abs(b') >= abs(crit95) & abs(bv')<crit99;
   Lcrit(1,:) = abs(b')>=crit99;
   
      
else;  % No significant pc-regression coefs -- null model
   clc, home;
   disp('NO SIGNIFICANT PC-REGRESSION COEFS IN MODEL');
	disp('ALL PCs TOGETHER EXPLAIN FOLLOWING PROPORTION TREE-RING VARIANCE: ');
   disp(V1b);
   disp('Press any key tO continue');
   pause;
   
end;  % of code dealing with Null model

%----- GATHER OUTPUT ARGUMENT DATA
datout{1}=V;
datout{2}=teach;
datout{3}=V1a;

%************** GRAPHICS CALL -- BAR GRAPH OF STANDARDIZED PCR WEIGHTS 
subfcn01(fwind,a',Lc',Clr,kcolor);

%****** GRAPHICS CALL-- BAR GRAPHS OF WEIGHTS ON ORIGINAL SEASONAL CLIMATE VARIABLES
fwind=fwind+1;

if nullfg2==0; % if at least one climate PC entered the PCR model
   subfcn03(fwind,dlab,slab,b',Lcrit,Clr,kcolor,strpd,treenm);
else; % No regression coefs significant in PCR model
   figure (fwind);
   title('NO RESPONSE FUNCTION JUSTIFIED: NO SIGNIFICANT PC PREDICTORS');
end




%********* SUBFUNCTIONS

function subfcn01(fwind,wgt,L,C,kopt)
%
% Sub-function of respfn1b.  PCR plots of stdzd weights
%
%************* IN
%
% fwind (1 x 1)i figure window
% wgt (1 x 2*nwgt)r  weights or correlations for the months or seas, P then T 
% L (3 x 2*nwgt)L  indicates whether weights signifant at .05, .10
%     row 1 applies to .01 level.  A "1" means signif, a "0" not
%     row 2 applies to .05 level. ....
%     row 3 applies to "insignificant" correls
% C (9 x 3)r  color schemes
%   Rows 1-3 for .01, .05 and other for color verion of P
%   Rows 4-6 for .01, ..... of T
%   Rows 7-9 for ....Black and White 
% kopt -----options
%   kopt(1)  color or b/w
%      ==1 color
%      ==2 b/w

figure(fwind)

statusc = close(fwind);
if statusc~=1;
   close(fwind);
end

figure(fwind);

set(gcf,'RendererMode','manual');
set(gcf,'Renderer','zbuffer');


set(gcf,'DefaultLineLineWidth',2.0);
set(gcf,'DefaultTextFontWeight','bold');

% Calculate number of weights
nwgt = length(wgt);
LP=L;

h1 = bar(wgt);
x1 = get(h1,'Xdata');
y1 = get(h1,'Ydata');


%Color patches
if kopt(1)==1; % color
   c1 =C(1,:); c2 = C(2,:); c3=C(3,:);
else
   c1 =C(7,:); c2=C(8,:); c3=C(9,:);
end
CP = zeros(1,nwgt,3);
n99 = sum(LP(1,:));
if n99>0;
   CP(1,LP(1,:),:) = repmat(c1,n99,1);
end
n95 = sum(LP(2,:));
if n95>0;
   CP(1,LP(2,:),:) = repmat(c2,n95,1);
end
nother = sum(LP(3,:));
if nother>0;
   CP(1,LP(3,:),:) = repmat(c3,nother,1);
end
hpatch1 = patch(x1,y1,CP);
set(hpatch1,'LineWidth',2);

% Adjust parent of patch
h1 = get(hpatch1,'parent');

set(h1,'XLim',[0 nwgt+2],...
   'XTick',[1:nwgt],...
   'Xgrid','on');

% Set + y limit at twice largest absolute value of any weight -- room for legend
maxwgt=max(wgt);
yylim = get(gca,'Ylim');
set(gca,'Ylim',[yylim(1)  2*maxwgt]);
yspace = maxwgt;

xlabel('Climate PC Number');


% Horiz Line at zero weight
hold on;
line([0 nwgt+1],[0 0],'Color',[0 0 0]);
hold off;


% Reference points for positioning legend
xxlim = get(gca,'XLim');
yylim = get(gca,'Ylim');
xwide = diff(xxlim);
ywide = diff(yylim);


%**************Color-patch legend
xgo = xwide/100;
xdel = xwide/30;
ygo1 = yylim(2)-yspace/20;
ydel = yspace/5;
ygo2 = ygo1 - 2*ydel;

% Upper part of legend
xleg =[xgo xgo  xgo+xdel xgo+xdel];
yleg = [ygo1-ydel ygo1 ygo1 ygo1-ydel];
hlegp1 = patch(xleg,yleg,c1);
text(xgo+xdel+xwide/100,ygo1-ydel/2,'\alpha=0.01','FontWeight','bold');

% Lower part of legend
xleg =[xgo xgo xgo+xdel xgo+xdel];
yleg = [ygo2-ydel ygo2 ygo2 ygo2-ydel];
hlegp2 = patch(xleg,yleg,c2);
text(xgo+xdel+xwide/100,ygo2-ydel/2,'\alpha=0.05');

ylabel('Standardized Coefficient');
title('PCR: Weights on PC''s of Seasonal Climate Variables',...
   'FontSize',12,'FontWeight','bold');



function subfcn03(fwind,dlab,slab,wgt,L,C,kopt,strpd,treenm)
%
% Sub-function for respfun1. PCR graphs of re-transformed weights on original
% seasonal climate variables
%
%************* IN
%
% fwind (1 x 1)i figure window
% dlab (2 x 1)s data labels (e.g., ['P';'T']
% slab (1 x ?) cell vector of season names
% wgt (1 x 2*nmos)r  weights on seasonalized P then T 
% L (3 x 2*nmos)L  indicates whether weights signifant at .05, .10
%     row 1 applies to .01 level.  A "1" means signif, a "0" not
%     row 2 applies to .05 level. ....
%     row 3 applies to "insignificant" correls
% C (9 x 3)r  color schemes
%   Rows 1-3 for .01, .05 and other for color verion of P
%   Rows 4-6 for .01, ..... of T
%   Rows 7-9 for ....Black and White 
% kopt -----options
%   kopt(1)  color or b/w
%      ==1 color
%      ==2 b/w
% strpd (1 x ?)s labeling string with period of analysis
% treenm (1 x ?s  label string of tree-ring series name
%
%*********** NOTES 
%
% 


figure(fwind)

statusc = close(fwind);
if statusc~=1;
   close(fwind);
end

figure(fwind);

set(gcf,'RendererMode','manual');
set(gcf,'Renderer','zbuffer');


set(gcf,'DefaultLineLineWidth',2.0);
set(gcf,'DefaultTextFontWeight','bold');

%Compute number of seasons
nmos = size(slab,2);


%----------  P axes (at bottom)
ax1 = axes;
set(ax1,'Position',[.1 .09 .8 .4]);

% pull P weights
pwgt = wgt (1:nmos);
LP = L(:,1:nmos);


h1 = bar(pwgt);
x1 = get(h1,'Xdata');
y1 = get(h1,'Ydata');


%Color patches
if kopt(1)==1; % color
   c1 =C(1,:); c2 = C(2,:); c3=C(3,:);
else
   c1 =C(7,:); c2=C(8,:); c3=C(9,:);
end
CP = zeros(1,nmos,3);
n99 = sum(LP(1,:));
if n99>0;
   CP(1,LP(1,:),:) = repmat(c1,n99,1);
end
n95 = sum(LP(2,:));
if n95>0;
   CP(1,LP(2,:),:) = repmat(c2,n95,1);
end
nother = sum(LP(3,:));
if nother>0;
   CP(1,LP(3,:),:) = repmat(c3,nother,1);
end
hpatch1 = patch(x1,y1,CP);
set(hpatch1,'LineWidth',2);

% Adjust parent of patch
h1 = get(hpatch1,'parent');

set(h1,'XLim',[0 nmos+2],...
   'XTick',[1:nmos],...
   'XTickLabel',slab,...
   'Xgrid','on');
%,...
  % 'YLim',[-1 1]);

xlabel('Season','FontWeight','bold');

xxlim = get(gca,'XLim');
xwide = diff(xxlim);
xgo = 0.93*xwide;

text(xgo,0.70,[dlab(1,:)],'FontSize',14);
set(gca,'FontWeight','bold');

hold on;
line([0 nmos+1],[0 0],'Color',[0 0 0]);
hold off;

%Color-patch legend
xgo = xwide/100;
xdel = xwide/30;
xleg =[xgo xgo  xgo+xdel xgo+xdel];
yleg = [.75 .95 .95 .75];
hlegp1 = patch(xleg,yleg,c1);
text(xgo+xdel+xwide/100,.85,'\alpha=0.01','FontWeight','bold');
xleg =[xgo xgo xgo+xdel xgo+xdel];
yleg = [.50 .70 .70 .50];
hlegp2 = patch(xleg,yleg,c2);
text(xgo+xdel+xwide/100,.60,'\alpha=0.05');

%************** T axes (one up from bottom)
ax2 = axes;
set(ax2,'Position',[.1 .52 .8 .4]);

% pull T weights
twgt = wgt ((nmos+1):(2*nmos));
LT = L(:,(nmos+1):(2*nmos));

hb = bar(twgt);
x2 = get(hb,'Xdata');
y2 = get(hb,'Ydata');

%Color patches
if kopt(1)==1; % color
   c1 =C(4,:); c2 = C(5,:); c3=C(6,:);
else
   c1 =C(7,:); c2=C(8,:); c3=C(9,:);
end
CT = zeros(1,nmos,3);
n99 = sum(LT(1,:));
if n99>0;
   CT(1,LT(1,:),:) = repmat(c1,n99,1);
end
n95 = sum(LT(2,:));
if n95>0;
   CT(1,LT(2,:),:) = repmat(c2,n95,1);
end
nother = sum(LT(3,:));
if nother>0;
   CT(1,LT(3,:),:) = repmat(c3,nother,1);
end

% Clear bar chart and add patch
cla;
hpatch2 = patch(x2,y2,CT);
set(hpatch2,'LineWidth',2);

% Adjust parent of patch
h2 = get(hpatch2,'parent');

set(h2,'XLim',[0 nmos+2],...
   'XTick',[1:nmos],...
   'XTickLabel',' ',...
   'Xgrid','on');

xxlim = get(gca,'XLim');
xwide = diff(xxlim);
xgo = 0.93*xwide;

text(xgo,0.75,[dlab(2,:)],'FontSize',14);


hold on;
line([0 nmos+1],[0 0],'color',[0 0 0]);
hold off;

set(gca,'FontWeight','bold');

%Color-patch legend
xgo = xwide/100;
xdel = xwide/30;
xleg =[xgo xgo  xgo+xdel xgo+xdel];
yleg = [.75 .95 .95 .75];
hlegt1 = patch(xleg,yleg,c1);
text(xgo+xdel+xwide/100,.85,'\alpha=0.01');
xleg =[xgo xgo xgo+xdel xgo+xdel];
yleg = [.50 .70 .70 .50];
hlegt2 = patch(xleg,yleg,c2);
text(xgo+xdel+xwide/100,.60,'\alpha=0.05');

txt1 = [treenm ';  ' strpd ';  PCR Weights on Seasonal Climate'];
title('PCR: Weights on Original Seasonal Climate Variables',...
   'FontSize',12,'FontWeight','bold');

