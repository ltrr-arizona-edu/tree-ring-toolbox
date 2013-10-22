function pbest=powtran(y,x,p,kmode)
% powtran: power transformation so that regression residuals not dependent on y
% pbest=powtran(y,x,p,kmode);
% Last revised 8-07-00
%
%*** IN
%
% y (my x 1)r dependent variable
% x (mx x 1)r independent variable, mx==my
% p (1 x ?)r  powers of x to try (e.g., p = (1:-0.5:0.2) in calling function)
% kmode (1 x 1)i rapid or informative mode

% Check sizes
[my,ny]=size(y);
[mx,nx]=size(x);
if my~=mx;
   error('x and y must be same size');
end;
if ny~=1 | nx~=1;
   error('x and y must be col vectors');
end;

% Check that x not negative
if any(x<0);
   error('all values in x must be non-negative');
end;

% Check that all powers greater than zero
if any (p<=0);
   error('p must be all positive');
end;


%*** BUILD MATRIX OF TRANSFORMED x

np =length(p); % number of powers
X=repmat(x,1,np); % dupe cv x into matrix
P=repmat(p,mx,1); % dupe rv p into matrix
T=X .^P;  % tsm of transformed x

%*** LOOP OVER POWERS OF X
Rsq1=repmat(NaN,np,1);
Rsq2=repmat(NaN,np,1);
for n = 1:np;
   str1=sprintf('%5.2f',p(n));
   % Regress y on x and compute Rsquared
   xthis=T(:,n);
   x1=[ones(mx,1) xthis];
   b=x1\y;
   yhat=x1 * b;
   e = y-yhat;
   Rsq = corrcoef(y,xthis);
   Rsq=Rsq(1,2) .^2; % prop var y explained
   Rsq1(n)=Rsq;
   str2=sprintf('%5.2f',Rsq);
   
   % Regress the residuals on y
   x1=[ones(mx,1) y]; 
   b=x1\e;
   ehat=x1 * b;
   Rsq = corrcoef(e,y);
   Rsq=Rsq(1,2) .^2; % prop var e explained
   Rsq2(n)=Rsq;
   str3=sprintf('%5.2f',Rsq);
   
   if kmode==2;
      subplot(1,2,1);
      plot(xthis,y,'o',xthis,yhat);
      title(['y on x for x raised to ' str1 '  R-sq = ' str2]);
      xlabel('Transformed x');
      ylabel('y');
      grid;
      %set(gca,'Position',[0.5780    0.1100    0.3270    0.8150]);
      subplot(1,2,2);
      plot(y,e,'o',y,ehat);
      title(['Residual on y.  Rsq = ' str3]);
      ylabel('e');
      xlabel('y');
      grid;
      %set(gca,'Position',[0.5780    0.1100    0.3270    0.8150]);
      pause;

   end;
end;
pbest=[];

      
        
   







